//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

#ifdef __CUDACC__
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

#ifdef __CUDACC__
  __device__ __constant__ short cHel[ncomb][npar];
  __device__ __constant__ fptype cIPC[%(ncouplingstimes2)i];
  __device__ __constant__ fptype cIPD[%(nparams)i];
  //__device__ __constant__ int cNGoodHel[1]; // FIXME: assume process.nprocesses == 1 for the moment
  __device__ __constant__ int cNGoodHel;
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
  static fptype cIPC[%(ncouplingstimes2)i];
  static fptype cIPD[%(nparams)i];
  //static int cNGoodHel[1]; // FIXME: assume process.nprocesses == 1 for the moment
  static int cNGoodHel;
  static int cGoodHel[ncomb];
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for given ihel
  // to running sum of |M|^2 over helicities for given event(s)
%(all_sigmaKin)s
  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( int numiterations,
                          int ngpublocks,
                          int ngputhreads,
                          bool verbose,
                          bool debug )
    : m_numiterations( numiterations )
    , gpu_nblocks( ngpublocks )
    , gpu_nthreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
    , m_pars( 0 )
    , m_masses()
    , dim( gpu_nblocks * gpu_nthreads )
  {
    // Helicities for the process - nodim
%(all_helicities)s
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(short) ) );
#else
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(short) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#ifndef __CUDACC__
    // SANITY CHECK: momenta AOSOA uses vectors with the same size as fptype_v
    //assert( neppV == mgOnGpu::neppM );
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
#ifdef __CUDACC__
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
    return out.str();
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            bool* isGoodHel )         // output: isGoodHel[ncomb] - device array
  {
    const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    fptype meHelSum[nprocesses] = { 0 }; // all zeros
    fptype meHelSumLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running
      // sum of |M|^2 over helicities for the given event
      calculate_wavefunctions( ihel, allmomenta, meHelSum[0] );
      if ( meHelSum[0] != meHelSumLast )
      {
        isGoodHel[ihel] = true;
        meHelSumLast = meHelSum[0];
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    //int nGoodHel[1] = { 0 }; // FIXME: assume process.nprocesses == 1 for the moment
    int nGoodHel = 0;
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
    //checkCuda( cudaMemcpyToSymbol( cNGoodHel, nGoodHel, sizeof(int) ) ); // FIXME: assume process.nprocesses == 1 for the moment
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof(int) ) );
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb * sizeof(int) ) );
  }
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2, part independent of incoming flavour
  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype_sv* allMEs            // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
  {
    mgDebugInitialise();

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

#ifndef __CUDACC__
    const int maxtry = 10;
    static unsigned long long sigmakin_itry = 0; // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = { false };
#endif
%(reset_jamp_lines)s
    // Start sigmaKin_lines
%(sigmaKin_lines)s
} // end namespace

//==========================================================================

// This was initially added to both C++ and CUDA in order to avoid RDC in CUDA (issue #51)
// This is now also needed by C++ LTO-like optimizations via inlining (issue #229)
%(hel_amps_cc)s

//==========================================================================
