//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 3.4.0_lo_vect, 2022-05-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <algorithm>
#include <iostream>
#include "mgOnGpuTypes.h"
#include "mgOnGpuConfig.h"

#include "CPPProcess.h"
#include "HelAmps_sm.h"
//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

//--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
template <typename mom_t, typename ipc_t, typename ipd_t>
KOKKOS_INLINE_FUNCTION fptype calculate_wavefunctions(
  const mom_t& allmomenta,              // input: momenta
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  fptype* allNumerators,                // output: multichannel numerators, running_sum_over_helicities
  fptype* allDenominators,              // output: multichannel denominators, running_sum_over_helicities
  const unsigned int channelId,         // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
  const short*  __restrict__ cHel,
  const ipc_t& COUPs,
  const ipd_t& cIPD
  )
{
  using namespace MG5_sm;
  fptype allMEs = 0;
  // The number of colors
  constexpr int ncolor = 1;

  // Local TEMPORARY variables for a subset of Feynman diagrams in the given event (ievt)
  // [NB these variables are reused several times (and re-initialised each time) within the same event or event page]
  cxtype w[mgOnGpu::nwf][mgOnGpu::nw6]; // particle wavefunctions within Feynman diagrams (nw6 is often 6, the dimension of spin 1/2 or spin 1 particles)
  cxtype amp[1]; // invariant amplitude for one given Feynman diagram

  // Local variables for the given event (ievt)
  cxtype jamp[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams in the event or event page

  // === Calculate wavefunctions and amplitudes for all diagrams in all processes - Loop over nevt events ===

  // Reset color flows (reset jamp) at the beginning of a new event or event page
  for( int i=0; i<ncolor; i++ ){ jamp[i] = cxmake( 0, 0 ); }

  // *** DIAGRAM 1 OF 2 ***

  // Wavefunction(s) for diagram number 1
#if not defined MGONGPU_TEST_DIVERGENCE
      opzxxx( Kokkos::subview(allmomenta, 0, Kokkos::ALL), cHel[0], -1, w[0] ); // NB: opzxxx only uses pz
#else
      if ( ievt % 2 == 0 )
        opzxxx( Kokkos::subview(allmomenta, 0, Kokkos::ALL), cHel[0], -1, w[0] ); // NB: opzxxx only uses pz
      else
        oxxxxx( Kokkos::subview(allmomenta, 0, 0, Kokkos::ALL), cHel[0], -1, w[0] );
#endif

  imzxxx( Kokkos::subview(allmomenta, 1, Kokkos::ALL), cHel[1], +1, w[1] ); // NB: imzxxx only uses pz

  ixzxxx( Kokkos::subview(allmomenta, 2, Kokkos::ALL), cHel[2], -1, w[2] );

  oxzxxx( Kokkos::subview(allmomenta, 3, Kokkos::ALL), cHel[3], +1, w[3] );

  FFV1P0_3( w[1], w[0], COUPs[0], 0., 0., w[4] );

  // Amplitude(s) for diagram number 1
  FFV1_0( w[2], w[3], w[4], COUPs[0], &amp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  if( channelId == 1 ) allNumerators[0] += cxabs2( amp[0] );
  if( channelId != 0 ) allDenominators[0] += cxabs2( amp[0] );
#endif
  jamp[0] -= amp[0];

  // *** DIAGRAM 2 OF 2 ***

  // Wavefunction(s) for diagram number 2
  FFV2_4_3( w[1], w[0], COUPs[1], COUPs[2], cIPD[0], cIPD[1], w[4] );

  // Amplitude(s) for diagram number 2
  FFV2_4_0( w[2], w[3], w[4], COUPs[1], COUPs[2], &amp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  if( channelId == 2 ) allNumerators[0] += cxabs2( amp[0] );
  if( channelId != 0 ) allDenominators[0] += cxabs2( amp[0] );
#endif
  jamp[0] -= amp[0];

  // *** COLOR ALGEBRA BELOW ***
  // (This method used to be called CPPProcess::matrix_1_epem_mupmum()?)

  // The color matrix

      // The color denominators (initialize all array elements, with ncolor=1)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype denom[ncolor] = { 1 }; // 1-D array[1]

      // The color matrix (initialize all array elements, with ncolor=1)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
      static constexpr fptype cf[ncolor][ncolor] = { { 1 } }; // 2-D array[1][1]

  // Sum and square the color flows to get the matrix element
  // (compute |M|^2 by squaring |M|, taking into account colours)
  fptype deltaMEs = { 0 }; // all zeros
  for( int icol = 0; icol < ncolor; icol++ )
  {
    cxtype ztemp = cxmake( 0, 0 );
    for( int jcol = 0; jcol < ncolor; jcol++ )
      ztemp += cf[icol][jcol] * jamp[jcol];
    deltaMEs += cxreal( ztemp * cxconj( jamp[icol] ) ) / denom[icol];
  }

  // *** STORE THE RESULTS ***

  // Store the leading color flows for choice of color
  // (NB: jamp2 must be an array of fptype)
  // for( int icol = 0; icol < ncolor; icol++ )
  // jamp2[0][icol] += cxreal( jamp[icol]*cxconj( jamp[icol] ) );

  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
  // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
  allMEs += deltaMEs;
  //printf( "calculate_wavefunction: %6d %2d %f\n", ievt, ihel, allMEs[ievt] );


  return allMEs;
}

//--------------------------------------------------------------------------

CPPProcess::CPPProcess( bool verbose, bool debug ):
        m_verbose(verbose),
        m_debug(debug)
        {
}

//--------------------------------------------------------------------------

// Initialize process (with parameters read from user cards)
void CPPProcess::initProc( const std::string& param_card_name ) {
    // Instantiate the model class and set parameters that stay fixed during run
    m_pars = Parameters_sm::getInstance();
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
    // Set external particle masses for this matrix element
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    m_masses.push_back( m_pars->ZERO );
    // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
    m_tIPC[0] = cxmake( m_pars->GC_3 );
    m_tIPC[1] = cxmake( m_pars->GC_50 );
    m_tIPC[2] = cxmake( m_pars->GC_59 );
    m_tIPD[0] = (fptype)m_pars->mdl_MZ;
    m_tIPD[1] = (fptype)m_pars->mdl_WZ;


}

// Define pointer accessors
cxtype* CPPProcess::get_tIPC_ptr() {return m_tIPC;}
const cxtype* CPPProcess::get_tIPC_ptr() const {return m_tIPC;}

fptype* CPPProcess::get_tIPD_ptr() {return m_tIPD;}
const fptype* CPPProcess::get_tIPD_ptr() const {return m_tIPD;}

bool kokkos_initialize() {
    //FIXME wrapper to Kokkos::initialize() here, might need modifications on certain platforms
    //// set Kokkos args
    //Kokkos::InitArguments args;
    ////args.device_id = 0;

    //// initialize Kokkos
    //Kokkos::initialize(args); 
    Kokkos::initialize();
    return true;
}

//--------------------------------------------------------------------------
template <typename mom_t, typename igh_t, typename ngh_t, 
          typename gs_t, typename idc_t, typename idp_t>
void sigmaKin_setup(
    const mom_t& momenta,
    const igh_t& iGoodHel,
    const ngh_t& nGoodHel,
    const gs_t& Gs,
    const idc_t& indep_couplings,
    const idp_t& cIPD,
    const int& league_size,
    const int& team_size)
{
  Kokkos::View<int*,Kokkos::DefaultExecutionSpace> isGoodHel("isGoodHel",mgOnGpu::ncomb); // has to be constant index, but should be `ncomb`

  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){
    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

    // Load helicities into local (private) memory
    //FIXME template variable does not work on NVIDIA devices
    //auto cHel = helicities<short>;
    constexpr short cHel[] {
      -1, -1, -1, -1,
      -1, -1, -1, 1,
      -1, -1, 1, -1,
      -1, -1, 1, 1,
      -1, 1, -1, -1,
      -1, 1, -1, 1,
      -1, 1, 1, -1,
      -1, 1, 1, 1,
      1, -1, -1, -1,
      1, -1, -1, 1,
      1, -1, 1, -1,
      1, -1, 1, 1,
      1, 1, -1, -1,
      1, 1, -1, 1,
      1, 1, 1, -1,
      1, 1, 1, 1
    };
    cxtype cIPC[dependentCouplings::ndcoup + independentCouplings::nicoup];

#if MGONGPU_NDCOUP > 0
    dependentCouplings::set_couplings_from_G(cIPC, Gs[ievt]); 
#endif

#if MGONGPU_NICOUP > 0
    for (size_t i = 0; i < independentCouplings::nicoup; i++) {
        cIPC[dependentCouplings::ndcoup + i] = indep_couplings(i);
    }
#endif


    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ihel = 0; ihel < mgOnGpu::ncomb; ihel++)
    {
      #ifdef MGONGPU_SUPPORTS_MULTICHANNEL
          constexpr unsigned int channelId = 0; // disable single-diagram channel enhancement
          fptype all_numerators = 0.0;
          fptype all_denominators = 0.0;
          auto allMEs = calculate_wavefunctions(local_mom, &all_numerators, &all_denominators, channelId, cHel + ihel*mgOnGpu::npar, cIPC, cIPD);
      #else
          auto allMEs = calculate_wavefunctions(local_mom, cHel + ihel*mgOnGpu::npar, cIPC, cIPD);
      #endif

      if (allMEs != 0)
      {
        isGoodHel(ihel) = true;
      }
    }
  });
  Kokkos::fence();

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,1),
  KOKKOS_LAMBDA(const int& i) {
      nGoodHel(0) = 0;
      for(int ihel=0; ihel < mgOnGpu::ncomb; ihel++) {
          if(isGoodHel(ihel)) {
              iGoodHel(nGoodHel(0)) = ihel;
              nGoodHel(0)++;
          }
      }
  });
  Kokkos::fence();
}


//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
template <typename mom_t, typename igh_t, typename ngh_t, 
          typename gs_t, typename idc_t, typename idp_t,
          typename out_t>
void sigmaKin(
    const mom_t& momenta,
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    const unsigned int channelId,          // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
    const igh_t& iGoodHel,
    const ngh_t& nGoodHel,
    const gs_t& Gs,
    const idc_t& indep_couplings,
    const idp_t& cIPD,
    const int& league_size,
    const int& team_size,
    const out_t& allMEs)
{
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    
    
    // Denominators: spins, colors and identical particles
    constexpr int denominators = 4; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Load helicities into local (private) memory
    //FIXME template variable does not work on NVIDIA devices
    //auto cHel = helicities<short>;
    constexpr short cHel[] {
      -1, -1, -1, -1,
      -1, -1, -1, 1,
      -1, -1, 1, -1,
      -1, -1, 1, 1,
      -1, 1, -1, -1,
      -1, 1, -1, 1,
      -1, 1, 1, -1,
      -1, 1, 1, 1,
      1, -1, -1, -1,
      1, -1, -1, 1,
      1, -1, 1, -1,
      1, -1, 1, 1,
      1, 1, -1, -1,
      1, 1, -1, 1,
      1, 1, 1, -1,
      1, 1, 1, 1
    };
    cxtype cIPC[dependentCouplings::ndcoup + independentCouplings::nicoup];

#if MGONGPU_NDCOUP > 0
    dependentCouplings::set_couplings_from_G(cIPC, Gs[ievt]); 
#endif

#if MGONGPU_NICOUP > 0
    for (size_t i = 0; i < independentCouplings::nicoup; i++) {
        cIPC[dependentCouplings::ndcoup + i] = indep_couplings(i);
    }
#endif

    allMEs[ievt] = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    fptype allNumerators = 0;
    fptype allDenominators = 0;
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (size_t ighel = 0; ighel < nGoodHel(0); ighel++) {
        const size_t ihel = iGoodHel[ighel];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      allMEs[ievt] += calculate_wavefunctions(local_mom, &allNumerators, &allDenominators, channelId, cHel + ihel*mgOnGpu::npar, cIPC, cIPD);
#else
      allMEs[ievt] += calculate_wavefunctions(local_mom, cHel + ihel*mgOnGpu::npar, cIPC, cIPD);
#endif
    }
    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) allMEs[ievt] *= allNumerators / allDenominators; // FIXME (#343): assume nprocesses == 1
#endif
    allMEs[ievt] /= (fptype)denominators;

  });// end parallel for
  Kokkos::fence();

}


//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


