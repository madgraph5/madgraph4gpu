//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

//--------------------------------------------------------------------------
%(all_sigmaKin)s
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
%(all_helicities)s
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
    constexpr int denominators = %(den_factors)s; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

    // Load helicities into local (private) memory
    //FIXME template variable does not work on NVIDIA devices
    //auto cHel = helicities<short>;
    constexpr short cHel[] {
%(all_helicities)s
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
