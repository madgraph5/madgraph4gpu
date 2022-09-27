//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

//--------------------------------------------------------------------------
%(all_sigmaKin)s
//--------------------------------------------------------------------------

CPPProcess::CPPProcess(
    int numiterations, int leaguesize, int teamsize,
    bool verbose, bool debug): 
      m_numiterations(numiterations), m_league_size(leaguesize), 
      m_team_size(teamsize),
      m_cHel("cHel",mgOnGpu::ncomb,mgOnGpu::npar), m_hHel("hHel",mgOnGpu::ncomb,mgOnGpu::npar), 
      m_cIPC("cIPC",%(ncouplingstimes2)s), m_hIPC("hIPC",%(ncouplingstimes2)s), 
      m_cIPD("cIPD",%(nparams)s), m_hIPD("hIPD",%(nparams)s)
{

  // Helicities for the process - nodim
  %(all_helicities)s

  for(int i=0;i<mgOnGpu::ncomb;++i)
      for(int j=0;j<mgOnGpu::npar;++j){
          m_hHel(i,j) = tHel[i][j];
      }
  Kokkos::deep_copy(m_cHel,m_hHel);
}

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
  // Read physics parameters like masses and couplings from user configuration files (static: initialize once)
%(assign_coupling)s

}

//--------------------------------------------------------------------------
template <typename mom_t, typename out_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin_setup(
    const mom_t& momenta,
    out_t& allMEs,
    const hel_t& cHel,
    const ipd_t& cIPD,
    const ipc_t& cIPC,
    igh_t& iGoodHel,
    ngh_t& nGoodHel,
    const int& ncomb,
    const int& league_size,
    const int& team_size) 
{
  Kokkos::View<int*,Kokkos::DefaultExecutionSpace> isGoodHel("isGoodHel",ncomb); // has to be constant index, but should be `ncomb`

  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){
    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    fptype allMEsLast = 0;

    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ihel = 0;ihel < ncomb;++ihel)
    {
      auto local_cHel = Kokkos::subview(cHel,ihel,Kokkos::ALL);
      calculate_wavefunctions(local_mom, local_cHel, cIPD, cIPC, allMEs[ievt]);
      if (allMEs[ievt] != allMEsLast)
      {
        isGoodHel(ihel) = true;
      }
      allMEsLast = allMEs[ievt];
    }
  });

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,1),
  KOKKOS_LAMBDA(const int& i){
    for(int ihel=0;ihel < ncomb;++ihel){

      if(isGoodHel(ihel)){
        iGoodHel(nGoodHel(0)) = ihel;
        nGoodHel(0)++;
      }

    }
  });
}


//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
template <typename mom_t, typename out_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin(const mom_t& momenta, out_t& allMEs, const hel_t& cHel,
    const ipd_t& cIPD, const ipc_t& cIPC, const igh_t& iGoodHel,
    const ngh_t& nGoodHel, const int& ncomb, const int& league_size,
    const int& team_size)
{
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    
    // Denominators: spins, colors and identical particles
    constexpr int denominators = 512;
    allMEs[ievt] = 0;

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ighel = 0;ighel < nGoodHel(0);++ighel)
    {
      auto local_cHel = Kokkos::subview(cHel,iGoodHel(ighel),Kokkos::ALL);
      calculate_wavefunctions(local_mom, local_cHel, cIPD, cIPC, allMEs[ievt]);
    }
    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    allMEs[ievt] /= (fptype)denominators;

  });// end parallel for

}


//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
