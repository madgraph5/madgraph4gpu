//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  using mgOnGpu::np4;   // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar;  // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
#ifdef MGONGPU_HARDCODE_CIPD
  %(hardcoded_assign_coupling)s
#else
#ifdef __CUDACC__
  __device__ __constant__ fptype cIPD[%(nparams)i];
#else
  static fptype cIPD[%(nparams)i];
#endif
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
#ifdef __CUDACC__
  __device__ __constant__ short cHel[ncomb][npar];
  __device__ __constant__ int cNGoodHel; // FIXME: assume process.nprocesses == 1 for the moment (eventually cNGoodHel[nprocesses]?)
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
  static int cNGoodHel; // FIXME: assume process.nprocesses == 1 for the moment (eventually cNGoodHel[nprocesses]?)
  static int cGoodHel[ncomb];
#endif

  //--------------------------------------------------------------------------
%(all_sigmaKin)s
  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( bool verbose,
                          bool debug )
    : m_verbose( verbose )
    , m_debug( debug )
#ifndef MGONGPU_HARDCODE_CIPD
    , m_pars( 0 )
#endif
    , m_masses()
  {
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
%(all_helicities)s
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * mgOnGpu::npar * sizeof( short ) ) );
#else
    memcpy( cHel, tHel, ncomb * mgOnGpu::npar * sizeof( short ) );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HARDCODE_CIPD
  // Initialize process (with parameters read from user cards)
  void
  CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_%(model_name)s::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    //m_pars->setDependentParameters(); // now computed event-by-event (running alphas #373)
    //m_pars->setDependentCouplings(); // now computed event-by-event (running alphas #373)
    if( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      //m_pars->printDependentParameters(); // now computed event-by-event (running alphas #373)
      //m_pars->printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    %(initProc_lines)s
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    %(assign_coupling)s
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, %(nparams)i * sizeof( fptype ) ) );
#else
    memcpy( cIPD, tIPD, %(nparams)i * sizeof( fptype ) );
#endif
    //for ( i=0; i<2; i++ ) std::cout << std::setprecision(17) << "tIPD[i] = " << tIPD[i] << std::endl;
  }
#else
  // Initialize process (with hardcoded parameters)
  void
  CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if( m_verbose )
    {
      Parameters_%(model_name)s::printIndependentParameters();
      Parameters_%(model_name)s::printIndependentCouplings();
      //Parameters_%(model_name)s::printDependentParameters(); // now computed event-by-event (running alphas #373)
      //Parameters_%(model_name)s::printDependentCouplings(); // now computed event-by-event (running alphas #373)
    }
    %(hardcoded_initProc_lines)s
  }
#endif

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string
  CPPProcess::getCompiler()
  {
    std::stringstream out;
    // CUDA version (NVCC)
    // [Use __NVCC__ instead of __CUDACC__ here!]
    // [This tests if 'nvcc' was used even to build a .cc file, even if not necessarily 'nvcc -x cu' for a .cu file]
    // [Check 'nvcc --compiler-options -dM -E dummy.c | grep CUDA': see https://stackoverflow.com/a/53713712]
#ifdef __NVCC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
    out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
    out << "nvcc UNKNOWN";
#endif
    out << " (";
#endif
    // ICX version (either as CXX or as host compiler inside NVCC)
#if defined __INTEL_COMPILER
#error "icc is no longer supported: please use icx"
#elif defined __INTEL_LLVM_COMPILER // alternative: __INTEL_CLANG_COMPILER
    out << "icx " << __INTEL_LLVM_COMPILER;
#ifdef __NVCC__
    out << ", ";
#else
    out << " (";
#endif
#endif
    // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
#ifdef __APPLE__
    out << "Apple clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#else
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
    // GCC toolchain version inside CLANG
    std::string tchainout;
    std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
    std::unique_ptr<FILE, decltype( &pclose )> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
    if( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
    std::array<char, 128> tchainbuf;
    while( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
    tchainout.pop_back(); // remove trailing newline
#if defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ", gcc " << tchainout;
#else
    out << " (gcc " << tchainout << ")";
#endif
#endif
#else
    out << "clang UNKNOWKN";
#endif
#else
    // GCC version (either as CXX or as host compiler inside NVCC)
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
    out << "gcc " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "gcc UNKNOWKN";
#endif
#endif
#if defined __NVCC__ or defined __INTEL_LLVM_COMPILER
    out << ")";
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

  __global__ void /* clang-format off */
  computeDependentCouplings( const fptype* allgs, // input: Gs[nevt]
                             fptype* allcouplings // output: couplings[nevt*ndcoup*2]
#ifndef __CUDACC__
                             , const int nevt     // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
  ) /* clang-format on */
  {
#ifdef __CUDACC__
    using namespace mg5amcGpu;
    using G_ACCESS = DeviceAccessGs;
    using C_ACCESS = DeviceAccessCouplings;
    G2COUP<G_ACCESS, C_ACCESS>( allgs, allcouplings );
#else
    using namespace mg5amcCpu;
    using G_ACCESS = HostAccessGs;
    using C_ACCESS = HostAccessCouplings;
    for( int ipagV = 0; ipagV < nevt / neppV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      const fptype* gs = MemoryAccessGs::ieventAccessRecordConst( allgs, ievt0 );
      fptype* couplings = MemoryAccessCouplings::ieventAccessRecord( allcouplings, ievt0 );
      G2COUP<G_ACCESS, C_ACCESS>( gs, couplings );
    }
#endif
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__ void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                       bool* isGoodHel )           // output: isGoodHel[ncomb] - device array
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs );
      if( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void
  sigmaKin_getGoodHel( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
                       const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
                       fptype* allMEs,             // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                       bool* isGoodHel,            // output: isGoodHel[ncomb] - device array
                       const int nevt )            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    //assert( (size_t)(allmomenta) %% mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) %% mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    const int maxtry0 = ( neppV > 16 ? neppV : 16 ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    fptype allMEsLast[maxtry0] = { 0 };              // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    const int maxtry = std::min( maxtry0, nevt );    // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for( int ievt = 0; ievt < maxtry; ++ievt )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, maxtry );
      for( int ievt = 0; ievt < maxtry; ++ievt )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        const bool differs = ( allMEs[ievt] != allMEsLast[ievt] );
        if( differs )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast[ievt] = allMEs[ievt]; // running sum up to helicity ihel
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  void
  sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel = 0;           // FIXME: assume process.nprocesses == 1 for the moment (eventually nGoodHel[nprocesses]?)
    int goodHel[ncomb] = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
    for( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if( isGoodHel[ihel] )
      {
        //goodHel[nGoodHel[0]] = ihel; // FIXME: assume process.nprocesses == 1 for the moment
        //nGoodHel[0]++; // FIXME: assume process.nprocesses == 1 for the moment
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof( int ) ) ); // FIXME: assume process.nprocesses == 1 for the moment
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb * sizeof( int ) ) );
#else
    cNGoodHel = nGoodHel;
    for( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  __global__ void /* clang-format off */
  sigmaKin( const fptype* allmomenta,   // input: momenta[nevt*npar*4]
            const fptype* allcouplings, // input: couplings[nevt*ndcoup*2]
            fptype* allMEs              // output: allMEs[nevt], |M|^2 final_avg_over_helicities
#ifndef __CUDACC__
            , const int nevt            // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
            ) /* clang-format on */
  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    constexpr int nprocesses = %(nproc)i;
    static_assert( nprocesses == 1, "Assume nprocesses == 1" ); // FIXME (#343): assume nprocesses == 1
    constexpr int denominators[%(nproc)i] = { %(den_factors)s };

#ifdef __CUDACC__
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
#else
    //assert( (size_t)(allmomenta) %% mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
    //assert( (size_t)(allMEs) %% mgOnGpu::cppAlign == 0 ); // SANITY CHECK: require SIMD-friendly alignment [COMMENT OUT TO TEST MISALIGNED ACCESS]
#endif

    // Start sigmaKin_lines
%(sigmaKin_lines)s
} // end namespace
