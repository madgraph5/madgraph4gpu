//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

namespace Proc
{
  static constexpr int np4 = mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  static constexpr int npar = mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  static constexpr int ncomb = mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  static constexpr int nwf = mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  static constexpr int nw6 = mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  static constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time

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
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  // Initialize process (with parameters read from user cards)
  void CPPProcess::initProc( const std::string& param_card_name )
  {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_%(model_name)s::getInstance();
    SLHAReader slha( param_card_name, m_verbose );
    m_pars->setIndependentParameters( slha );
    m_pars->setIndependentCouplings();
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();
    if ( m_verbose )
    {
      m_pars->printIndependentParameters();
      m_pars->printIndependentCouplings();
      //m_pars->printDependentParameters();
      //m_pars->printDependentCouplings();
    }
    %(initProc_lines)s
#ifndef MGONGPU_HARDCODE_PARAM
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
%(assign_coupling)s
#endif
  }

  //--------------------------------------------------------------------------
#ifndef MGONGPU_HARDCODE_PARAM
  // Define pointer accessors
  cxtype* CPPProcess::get_tIPC_ptr() {return m_tIPC;}
  const cxtype* CPPProcess::get_tIPC_ptr() const {return m_tIPC;}

  fptype* CPPProcess::get_tIPD_ptr() {return m_tIPD;}
  const fptype* CPPProcess::get_tIPD_ptr() const {return m_tIPD;}
#endif

  //--------------------------------------------------------------------------

  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                            bool* isGoodHel,                       // output: isGoodHel[ncomb] - device array
                            const short* __restrict__ cHel,
                            const cxtype* __restrict__ cIPC,
                            const fptype* __restrict__ cIPD
                            )
  {
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    fptype allMEs = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      constexpr unsigned int channelId = 0; // disable single-diagram channel enhancement
      fptype allNumerators = 0;
      fptype allDenominators = 0;
      allMEs += calculate_wavefunctions( allmomenta, &allNumerators, &allDenominators, channelId, cHel + ihel*npar, cIPC, cIPD );
#else
      allMEs += calculate_wavefunctions( allmomenta, cHel + ihel*npar, cIPC, cIPD );
#endif
      if ( allMEs != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs; // running sum up to helicity ihel for event ievt
    }
  }

  //--------------------------------------------------------------------------

  int sigmaKin_setGoodHel( const bool* isGoodHel, int* goodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel = 0; // FIXME: assume process.nprocesses == 1 for the moment (eventually nGoodHel[nprocesses]?)
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
    return nGoodHel;
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  SYCL_EXTERNAL
  fptype sigmaKin( const fptype* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                   const unsigned int channelId,          // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
                   const short* __restrict__ cHel,
                   const cxtype* __restrict__ cIPC,
                   const fptype* __restrict__ cIPD,
                   const int* __restrict__ cNGoodHel,
                   const int* __restrict__ cGoodHel
                 )
  {
    mgDebugInitialise();

    // Denominators: spins, colors and identical particles
    constexpr int denominators = %(den_factors)s; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    //m_pars->setDependentParameters();
    //m_pars->setDependentCouplings();

    // Start sigmaKin_lines
%(sigmaKin_lines)s
} // end namespace

//==========================================================================
