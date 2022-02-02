//==========================================================================
// Class member functions for calculating the matrix elements for
%(process_lines)s

using mgKokkos::ncomb; // number of helicity combinations
//--------------------------------------------------------------------------

// Evaluate |M|^2 for each subprocess
// NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event


%(all_sigmaKin)s

//--------------------------------------------------------------------------

template <typename mom_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin_setup(
    const mom_t& momenta,
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
  KOKKOS_LAMBDA(member_type team_member){
    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

    const int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
    double meHelSum[nprocesses] = {0};  // all zeros
    double meHelSumLast = 0;

    auto evt_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ihel = 0;ihel < ncomb;++ihel)
    {
      auto local_cHel = Kokkos::subview(cHel,ihel,Kokkos::ALL);
      calculate_wavefunctions(local_cHel, evt_mom, cIPD, cIPC, meHelSum[0]);
      if (meHelSum[0] != meHelSumLast)
      {
        isGoodHel(ihel) = true;
        meHelSumLast = meHelSum[0];
      }
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
  KOKKOS_LAMBDA(member_type team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    
    // Denominators: spins, colors and identical particles
    constexpr int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
    constexpr int denominators[1] = {512};

    // Reset the "matrix elements" - running sums of |M|^2 over helicities for
    // the given event
    double meHelSum[nprocesses] = {0};  // all zeros

    // using precomputed good helicities
    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ighel = 0;ighel < nGoodHel(0);++ighel)
    {
      auto local_cHel = Kokkos::subview(cHel,iGoodHel(ighel),Kokkos::ALL);
      calculate_wavefunctions(local_cHel, local_mom, cIPD, cIPC, meHelSum[0]);
    }
    // Get the final |M|^2 as an average over helicities/colors of the running
    // sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf
    // 1d7c/Handout_4_2016-UZH.pdf]
    for (int iproc = 0;iproc < nprocesses; ++iproc)
    {
      meHelSum[iproc] /= denominators[iproc];
    }

    // Set the final average |M|^2 for this event in the output array for all
    // events
    for (int iproc = 0;iproc < nprocesses; ++iproc)
    {
      allMEs[iproc * nprocesses + ievt] = meHelSum[iproc];
    }
  });// end parallel for

}


//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
