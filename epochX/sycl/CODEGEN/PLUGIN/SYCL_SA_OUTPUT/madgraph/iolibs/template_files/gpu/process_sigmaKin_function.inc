    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEs = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    fptype allNumerators = 0;
    fptype allDenominators = 0;
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    for ( int ighel = 0; ighel < cNGoodHel[0]; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      allMEs += calculate_wavefunctions( allmomenta, &allNumerators, &allDenominators, channelId, cHel + ihel*npar, cIPC, cIPD );
#else
      allMEs += calculate_wavefunctions( allmomenta, cHel + ihel*npar, cIPC, cIPD );
#endif
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    mgDebugFinalise();
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) allMEs *= allNumerators / allDenominators; // FIXME (#343): assume nprocesses == 1
#endif
    return allMEs / denominators;
  }

  //--------------------------------------------------------------------------
