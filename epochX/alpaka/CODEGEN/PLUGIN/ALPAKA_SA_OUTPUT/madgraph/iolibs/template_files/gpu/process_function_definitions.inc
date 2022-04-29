//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

#if defined(ALPAKA) || defined(__CUDACC__)
namespace gProc
#else
namespace Proc
#endif
{
  using mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  // Physics parameters (masses, coupling, etc...)
  // For CUDA performance, hardcoded constexpr's would be better: fewer registers and a tiny throughput increase
  // However, physics parameters are user-defined through card files: use CUDA constant memory instead (issue #39)
  // [NB if hardcoded parameters are used, it's better to define them here to avoid silent shadowing (issue #263)]
  //constexpr fptype cIPC[%(ncouplingstimes2)i] = { ... };
  //constexpr fptype cIPD[%(nparams)i] = { ... };
#ifdef ALPAKA
  ALPAKA_STATIC_ACC_MEM_CONSTANT fptype cIPC[%(ncouplingstimes2)i];
  ALPAKA_STATIC_ACC_MEM_CONSTANT fptype cIPD[%(nparams)i];
#elif defined(__CUDACC__)
  __device__ __constant__ fptype cIPC[%(ncouplingstimes2)i];
  __device__ __constant__ fptype cIPD[%(nparams)i];
#else
  static fptype cIPC[%(ncouplingstimes2)i];
  static fptype cIPD[%(nparams)i];
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
#ifdef ALPAKA
  ALPAKA_STATIC_ACC_MEM_CONSTANT short cHel[ncomb][npar];
  ALPAKA_STATIC_ACC_MEM_CONSTANT int cNGoodHel;
  ALPAKA_STATIC_ACC_MEM_CONSTANT int cGoodHel[ncomb];
#elif defined(__CUDACC__)
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

  CPPProcess::CPPProcess( int numiterations,
                          int ngpublocks,
                          int ngputhreads,
                          bool verbose,
                          bool debug )
    : m_numiterations( numiterations )
    , m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
    , m_pars( 0 )
    , m_masses()
  {
    // Helicities for the process - nodim
%(all_helicities)s
#ifdef ALPAKA
    auto extent = alpaka::Vec<alpaka::DimInt<1u>, int>{ncomb * mgOnGpu::npar};
    auto& device(cupla::manager::Device<cupla::AccDev>::get().current());
    auto viewcHel = alpaka::createStaticDevMemView((short*)cHel, device, extent);
    auto& streamObject(cupla::manager::Stream<cupla::AccDev,cupla::AccStream>::get().stream( 0 ));
    auto& host(cupla::manager::Device<cupla::AccHost>::get().current());
    alpaka::ViewPlainPtr<cupla::AccHost, short, alpaka::DimInt<1u>, int> bufHostHel((short*)tHel, host, extent);
    alpaka::memcpy(streamObject, viewcHel, bufHostHel, extent);
    alpaka::wait( streamObject );
#elif defined(__CUDACC__)
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * mgOnGpu::npar * sizeof(short) ) );
#else
    memcpy( cHel, tHel, ncomb * mgOnGpu::npar * sizeof(short) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#if !defined(ALPAKA) && !defined(__CUDACC__)
    // SANITY CHECK: momenta AOSOA uses vectors with the same size as fptype_v
    assert( neppV == mgOnGpu::neppM );
#endif
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  // Initialize process
  void CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    if ( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
    }
    m_pars->setDependentParameters();
    m_pars->setDependentCouplings();

    %(initProc_lines)s

    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    %(assign_coupling)s
#ifdef ALPAKA
    auto extentcx3 = alpaka::Vec<alpaka::DimInt<1u>, int>{6};
    auto extentfp2 = alpaka::Vec<alpaka::DimInt<1u>, int>{2};
    auto& device(cupla::manager::Device<cupla::AccDev>::get().current());
    auto viewcIPC = alpaka::createStaticDevMemView((fptype*)cIPC, device, extentcx3);
    auto viewcIPD = alpaka::createStaticDevMemView((fptype*)cIPD, device, extentfp2);
    auto& streamObject(cupla::manager::Stream<cupla::AccDev,cupla::AccStream>::get().stream( 0 ));
    auto& host(cupla::manager::Device<cupla::AccHost>::get().current());
    alpaka::ViewPlainPtr<cupla::AccHost, fptype, alpaka::DimInt<1u>, int> bufHostIPC((fptype *)tIPC, host, extentcx3);
    alpaka::ViewPlainPtr<cupla::AccHost, fptype, alpaka::DimInt<1u>, int> bufHostIPD((fptype *)tIPD, host, extentfp2);
    alpaka::memcpy(streamObject, viewcIPC, bufHostIPC, extentcx3);
    alpaka::wait( streamObject );
    alpaka::memcpy(streamObject, viewcIPD, bufHostIPD, extentfp2);
    alpaka::wait( streamObject );
#elif defined(__CUDACC__)
    checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, %(ncouplings)i * sizeof(cxtype) ) );
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, %(nparams)i * sizeof(fptype) ) );
#else
    memcpy( cIPC, tIPC, %(ncouplings)i * sizeof(cxtype) );
    memcpy( cIPD, tIPD, %(nparams)i * sizeof(fptype) );
#endif

    //std::cout << std::setprecision(17) << "tIPC[0] = " << tIPC[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[1] = " << tIPC[1] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[2] = " << tIPC[2] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[0] = " << tIPD[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[1] = " << tIPD[1] << std::endl;
  }

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string CPPProcess::getCompiler()
  {
    std::stringstream out;
    // CUDA version (NVCC)
#ifdef __CUDACC__
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
    out << "icx " << __INTEL_LLVM_COMPILER << " (";
#endif
    // CLANG version (either as CXX or as host compiler inside NVCC or inside ICX)
#if defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
    // GCC toolchain version inside CLANG
    std::string tchainout;
    std::string tchaincmd = "readelf -p .comment $(${CXX} -print-libgcc-file-name) |& grep 'GCC: (GNU)' | grep -v Warning | sort -u | awk '{print $5}'";
    std::unique_ptr<FILE, decltype(&pclose)> tchainpipe( popen( tchaincmd.c_str(), "r" ), pclose );
    if ( !tchainpipe ) throw std::runtime_error( "`readelf ...` failed?" );
    std::array<char, 128> tchainbuf;
    while ( fgets( tchainbuf.data(), tchainbuf.size(), tchainpipe.get() ) != nullptr ) tchainout += tchainbuf.data();
    tchainout.pop_back(); // remove trailing newline
#if defined __CUDACC__ or defined __INTEL_LLVM_COMPILER
    out << ", gcc " << tchainout;
#else
    out << " (gcc " << tchainout << ")";
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
#if defined __CUDACC__ or defined __INTEL_LLVM_COMPILER
    out << ")";
#endif
#ifdef ALPAKA
    out << " (ALPAKA)";
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

#if defined(ALPAKA) || defined(__CUDACC__)
#ifdef ALPAKA
template< typename T_Acc >
  ALPAKA_FN_ACC
  void sigmaKin_getGoodHel::operator()( T_Acc const &acc,
                            const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel ) const      // output: isGoodHel[ncomb] - device array
#else
  __global__
  void sigmaKin_getGoodHel( 
                            const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel )            // output: isGoodHel[ncomb] - device array
#endif
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
#ifdef ALPAKA
      calculate_wavefunctions( acc, ihel, allmomenta, allMEs );
#else
      calculate_wavefunctions( ihel, allmomenta, allMEs );
#endif
      if ( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype_sv* allMEs,           // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                            bool* isGoodHel              // output: isGoodHel[ncomb] - device array
                            , const int nevt )           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int maxtry0 = ( neppV > 16 ? neppV : 16 ); // 16, but at least neppV (otherwise the npagV loop does not even start)
    fptype_sv allMEsLast[maxtry0/neppV] = { 0 };
    const int maxtry = std::min( maxtry0, nevt ); // 16, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for ( int ipagV = 0; ipagV < maxtry/neppV; ++ipagV )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ipagV] = fptype_sv{0}; // all zeros
    }
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      calculate_wavefunctions( ihel, allmomenta, allMEs, maxtry );
      for ( int ipagV = 0; ipagV < maxtry/neppV; ++ipagV )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        const bool differs = maskor( allMEs[ipagV] != allMEsLast[ipagV] ); // true if any of the neppV events differs
        if ( differs )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast[ipagV] = allMEs[ipagV]; // running sum up to helicity ihel
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel = 0; // FIXME: assume process.nprocesses == 1 for the moment (eventually nGoodHel[nprocesses]?)
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        //goodHel[nGoodHel[0]] = ihel; // FIXME: assume process.nprocesses == 1 for the moment
        //nGoodHel[0]++; // FIXME: assume process.nprocesses == 1 for the moment
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef ALPAKA
    auto extent1 = alpaka::Vec<alpaka::DimInt<1u>, int>{1};
    auto extentn = alpaka::Vec<alpaka::DimInt<1u>, int>{ncomb};
    auto& device(cupla::manager::Device<cupla::AccDev>::get().current());
    auto viewcNgood = alpaka::createStaticDevMemView(&cNGoodHel, device, extent1);
    auto viewcGood = alpaka::createStaticDevMemView((int*)cGoodHel, device, extentn);
    auto& streamObject(cupla::manager::Stream<cupla::AccDev,cupla::AccStream>::get().stream( 0 ));
    auto& host(cupla::manager::Device<cupla::AccHost>::get().current());
    alpaka::ViewPlainPtr<cupla::AccHost, int, alpaka::DimInt<1u>, int> bufHostNgood(&nGoodHel, host, extent1);
    alpaka::ViewPlainPtr<cupla::AccHost, int, alpaka::DimInt<1u>, int> bufHostGood((int*)goodHel, host, extentn);
    alpaka::memcpy(streamObject, viewcNgood, bufHostNgood, extent1);
    alpaka::wait( streamObject );
    alpaka::memcpy(streamObject, viewcGood, bufHostGood, extentn);
    alpaka::wait( streamObject );
#elif defined(__CUDACC__)
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof(int) ) ); // FIXME: assume process.nprocesses == 1 for the moment
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb*sizeof(int) ) );
#else
    cNGoodHel = nGoodHel;
    for ( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

#ifdef ALPAKA
  template< typename T_Acc >
  ALPAKA_FN_ACC
  void sigmaKin::operator()( T_Acc const &acc,
                             const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                             fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
                           ) const
#else
  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype_sv* allMEs            // output: allMEs[npagM][neppM], final |M|^2 averaged over helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
#endif
  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    const int denominators = %(den_factors)s; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

#if defined(ALPAKA) || defined(__CUDACC__)
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "sigmakin: ievt %%d\n", ievt );
#endif

    // Start sigmaKin_lines
%(sigmaKin_lines)s
} // end namespace

//==========================================================================

// This was initially added to both C++ and CUDA in order to avoid RDC in CUDA (issue #51)
// This is now also needed by C++ LTO-like optimizations via inlining (issue #229)
%(hel_amps_cc)s

//==========================================================================

#ifdef ALPAKA
template ALPAKA_FN_ACC void gProc::sigmaKin_getGoodHel::operator()<cupla::Acc>(cupla::Acc const&, const fptype_sv* , fptype_sv*, bool *) const;
template ALPAKA_FN_ACC void gProc::sigmaKin::operator()<cupla::Acc>(cupla::Acc const&, const fptype_sv* , fptype_sv*) const;
#endif
