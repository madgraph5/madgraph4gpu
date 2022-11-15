 
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(const member_type& team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)

    auto local_me = Kokkos::subview(allMEs,ievt);
    local_me() = 0;
    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    auto local_numerator = Kokkos::subview(allNumerators,ievt);
    auto local_denomator = Kokkos::subview(allDenominators,ievt);
    local_numerator() = 0;
    local_denomator() = 0;
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    for (int ighel = 0; ighel < nGoodHel(0); ighel++)
    {
      auto local_cHel = Kokkos::subview(cHel,iGoodHel(ighel),Kokkos::ALL);
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      local_me() += calculate_wavefunctions(local_mom, local_numerator(), local_denomator(), channelId, local_cHel, cIPC, cIPD);
#else
      local_me() += calculate_wavefunctions(local_mom, local_cHel, cIPC, cIPD);
#endif
    }
    
    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    mgDebugFinalise();
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) local_me() *= local_numerator() / local_denomator(); // FIXME (#343): assume nprocesses == 1
#endif
    local_me() /= local_denomator();

  });// end parallel for
  //--------------------------------------------------------------------------
