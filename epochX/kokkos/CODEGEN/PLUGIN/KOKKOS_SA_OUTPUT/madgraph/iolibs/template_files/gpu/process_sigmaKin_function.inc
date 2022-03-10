   
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(member_type team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    
    // Denominators: spins, colors and identical particles
    constexpr int nprocesses = %(nproc)i;  // FIXME: assume process.nprocesses == 1
    constexpr int denominators[%(nproc)i] = {%(den_factors)s};

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

 
