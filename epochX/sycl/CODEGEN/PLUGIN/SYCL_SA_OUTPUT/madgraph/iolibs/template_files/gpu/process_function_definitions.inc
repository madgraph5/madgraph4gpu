//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

namespace Proc
{
  using mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  using mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  using mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)

  // Physics parameters (masses, coupling, etc...)
#ifdef MGONGPU_HARDCODE_CIPC
  %(assign_hardcoded_coupling)s
#else
  static fptype cIPC[%(ncouplingstimes2)i];
  static fptype cIPD[%(nparams)i];
#endif

  // Helicity combinations (and filtering of "good" helicity combinations)
  static int cHel[ncomb][npar];

  //--------------------------------------------------------------------------
%(all_sigmaKin)s
  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess(
          int numiterations,
          int ngpublocks,
          int ngputhreads,
          bool verbose,
          bool debug )
    : m_numiterations( numiterations )
    , m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
#ifndef MGONGPU_HARDCODE_CIPC
    , m_pars( 0 )
#endif
    , m_masses()
  {
    // Helicities for the process [NB do keep 'static' for this constexpr array, see issue #283]
%(all_helicities)s
    memcpy( *cHel, *tHel, ncomb * mgOnGpu::npar * sizeof(short) );

    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    static_assert( sizeof(cxtype) == 2 * sizeof(fptype) );
    // SANITY CHECK: check that neppR, neppM and neppV are powers of two (https://stackoverflow.com/a/108360)
    auto ispoweroftwo = []( int n ) { return ( n > 0 ) && !( n & ( n - 1 ) ); };
    static_assert( ispoweroftwo( mgOnGpu::neppR ) );
    static_assert( ispoweroftwo( mgOnGpu::neppM ) );
    static_assert( ispoweroftwo( neppV ) );
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  //const std::vector<fptype>& CPPProcess::getMasses() const { return mME; }

  //--------------------------------------------------------------------------

  short * CPPProcess::get_cHel_ptr() const {return *cHel;}

  //--------------------------------------------------------------------------

  fptype* CPPProcess::get_cIPC_ptr() {return cIPC;}

  //--------------------------------------------------------------------------

  fptype* CPPProcess::get_cIPD_ptr() {return cIPD;}

  //--------------------------------------------------------------------------

#ifndef MGONGPU_HARDCODE_CIPC
  // Initialize process (with parameters read from user cards)
  void CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    m_pars->setDependentParameters();
    m_pars->setDependentCouplings();
    if ( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      m_pars->printDependentParameters();
      m_pars->printDependentCouplings();
    }
    %(initProc_lines)s
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    // Then copy them to CUDA constant memory (issue #39) or its C++ emulation in file-scope static memory
    %(assign_coupling)s
    memcpy( cIPC, tIPC, %(ncouplings)i * sizeof(cxtype) );
    memcpy( cIPD, tIPD, %(nparams)i * sizeof(fptype) );
  }
#else
  // Initialize process (with hardcoded parameters)
  void CPPProcess::initProc( const std::string& /*param_card_name*/ )
  {
    // Use hardcoded physics parameters
    if ( m_verbose )
    {
      Parameters_sm::printIndependentParameters();
      Parameters_sm::printIndependentCouplings();
      Parameters_sm::printDependentParameters();
      Parameters_sm::printDependentCouplings();
    }
    // Set external particle masses for this matrix element
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
    m_masses.push_back( Parameters_sm::ZERO );
  }
#endif

  //--------------------------------------------------------------------------

  // Retrieve the compiler that was used to build this module
  const std::string CPPProcess::getCompiler() {
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

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,           // output: allMEs[nevt], |M|^2 final_avg_over_helicities
                            bool * isGoodHel,  // output: isGoodHel[ncomb] - device array
                            size_t ievt,
                            int *cHel,
                            const fptype *cIPC,
                            const fptype *cIPD
                            ) {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype allMEsLast = 0;
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
        calculate_wavefunctions( ihel, allmomenta, allMEs, ievt, cHel, cIPC, cIPD );
        if ( allMEs[ievt] != allMEsLast )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
      }
  }

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_setGoodHel( const bool* isGoodHel, int* cNGoodHel, int* cGoodHel) {
    *cNGoodHel = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      if ( isGoodHel[ihel] ) {
        cGoodHel[*cNGoodHel] = ihel;
        (*cNGoodHel)++;
      }
    }
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  SYCL_EXTERNAL
  void sigmaKin(
          const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
          fptype* allMEs,            // output: allMEs[nevt], |M|^2 final_avg_over_helicities
          size_t ievt,
          int *cHel,
          const fptype *cIPC,
          const fptype *cIPD,
          int *cNGoodHel,
          int *cGoodHel
          )  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    const int denominators = %(den_factors)s; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

    // Start sigmaKin_lines
%(sigmaKin_lines)s
} // end namespace

