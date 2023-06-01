! Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
! Created by: J. Alwall (Jul 2010) for the MG5aMC CPP backend.
!==========================================================================
! Copyright (C) 2020-2023 CERN and UCLouvain.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
! Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
!==========================================================================
! Copyright (C) 2021-2023 Argonne National Laboratory.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
!==========================================================================
//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

namespace Proc
{
  static constexpr size_t np4 = mgOnGpu::np4; // dimensions of 4-momenta (E,px,py,pz)
  static constexpr size_t npar = mgOnGpu::npar; // #particles in total (external = initial + final): e.g. 4 for e+ e- -> mu+ mu-
  static constexpr size_t ncomb = mgOnGpu::ncomb; // #helicity combinations: e.g. 16 for e+ e- -> mu+ mu- (2**4 = fermion spin up/down ** npar)

  static constexpr size_t nwf = mgOnGpu::nwf; // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  static constexpr size_t nw6 = mgOnGpu::nw6; // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  static constexpr size_t neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time

  // The number of colors
  static constexpr size_t ncolor = %(ncolor)s;

  //--------------------------------------------------------------------------
%(all_sigmaKin)s
  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( size_t numiterations,
                          size_t ngpublocks,
                          size_t ngputhreads,
                          bool verbose,
                          bool debug )
    : m_numiterations( numiterations )
    , m_ngpublocks( ngpublocks )
    , m_ngputhreads( ngputhreads )
    , m_verbose( verbose )
    , m_debug( debug )
    #ifndef MGONGPU_HARDCODE_PARAM
        , m_pars( 0 )
        , m_masses()
    #endif
  {
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  // Initialize process (with parameters read from user cards)
  void CPPProcess::initProc( const std::string& param_card_name ) {
      #ifndef MGONGPU_HARDCODE_PARAM
          // Instantiate the model class and set parameters that stay fixed during run
          m_pars = Parameters_%(model_name)s::getInstance();
          SLHAReader slha( param_card_name, m_verbose );
          m_pars->setIndependentParameters( slha );
          m_pars->setIndependentCouplings();
          //m_pars->setDependentParameters();
          //m_pars->setDependentCouplings();
          if ( m_verbose ) {
              m_pars->printIndependentParameters();
              m_pars->printIndependentCouplings();
              //m_pars->printDependentParameters();
              //m_pars->printDependentCouplings();
          }
          %(initProc_lines)s
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
  void sigmaKin_getGoodHel( const vector4* __restrict__ allmomenta, // input: momenta[nevt*npar*4]
                            bool* isGoodHel,                        // output: isGoodHel[ncomb] - device array
                            const signed char* __restrict__ cHel,
                            const cxtype_sv* __restrict__ COUPs,
                            const fptype* __restrict__ cIPD
                            ) {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype_sv allMEsLast = FPZERO_SV;
      fptype_sv allMEs = FPZERO_SV;
      for ( size_t ihel = 0; ihel < ncomb; ihel++ ) {
          // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
          constexpr fptype_sv* jamp2_sv = nullptr; // no need for color selection during helicity filtering
          #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
              constexpr size_t channelId = 0; // disable single-diagram channel enhancement
              fptype_sv allNumerators = FPZERO_SV;
              fptype_sv allDenominators = FPZERO_SV;
              allMEs += calculate_wavefunctions( allmomenta, &allNumerators, &allDenominators, channelId, cHel + ihel*npar, COUPs, cIPD, jamp2_sv );
          #else
              allMEs += calculate_wavefunctions( allmomenta, cHel + ihel*npar, COUPs, cIPD, jamp2_sv );
          #endif
          if (FPANY_SV(allMEs != allMEsLast)) {
              isGoodHel[ihel] = true;
          }
          allMEsLast = allMEs; // running sum up to helicity ihel for event ievt
       }
   }

  //--------------------------------------------------------------------------

  size_t sigmaKin_setGoodHel( const bool* isGoodHel, size_t* goodHel ) {
      size_t nGoodHel = 0; // FIXME: assume process.nprocesses == 1 for the moment (eventually nGoodHel[nprocesses]?)
      for (size_t ihel = 0; ihel < ncomb; ihel++) {
          if (isGoodHel[ihel]) {
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
  fptype_sv sigmaKin( const vector4* __restrict__ allmomenta, // input: momenta[]
                      const fptype_sv* __restrict__ rndhel,   // input: random numbers[] for helicity selection
                      const fptype_sv* __restrict__ rndcol,   // input: random numbers[] for color selection
                      int_sv* __restrict__ selhel,            // output: helicity selection[]
                      int_sv* __restrict__ selcol,            // output: color selection[]
                      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                          const size_t channelId,             // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
                      #endif
                      const signed char* __restrict__ cHel,
                      const cxtype_sv* __restrict__ COUPs,
                      const fptype* __restrict__ cIPD,
                      const size_t* __restrict__ cNGoodHel,
                      const size_t* __restrict__ cGoodHel
                 ) {

      // Denominators: spins, colors and identical particles
      constexpr int denominators = %(den_factors)s; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

      // Start sigmaKin_lines
%(sigmaKin_lines)s
} // end namespace

//==========================================================================
