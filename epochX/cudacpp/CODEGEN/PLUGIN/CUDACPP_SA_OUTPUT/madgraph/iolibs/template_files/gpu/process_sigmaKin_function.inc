#ifndef __CUDACC__
    const int maxtry = 10;
    static unsigned long long sigmakin_itry = 0; // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = { false };
#endif

#ifndef __CUDACC__
    for ( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype_sv& allMEsThisEvt = allMEs[ievt];
      allMEsThisEvt = 0;

#ifdef __CUDACC__
      // CUDA - using precomputed good helicities
      for ( int ighel = 0; ighel < cNGoodHel; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        calculate_wavefunctions( ihel, allmomenta, &allMEsThisEvt );
      }
#else
      // C++ - compute good helicities within this loop
      fptype allMEsThisEvtLast = 0; // check for good helicities
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry > maxtry && !sigmakin_goodhel[ihel] ) continue;
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
        calculate_wavefunctions( ihel, allmomenta, &allMEsThisEvt, ievt );
        if ( sigmakin_itry <= maxtry )
        {
          if ( !sigmakin_goodhel[ihel] && allMEsThisEvt > allMEsThisEvtLast ) sigmakin_goodhel[ihel] = true;
          allMEsThisEvtLast = allMEsThisEvt;
        }
      }
#endif

      // Set the final average |M|^2 for this event in the output array for all events
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = allMEsThisEvt;

#ifndef __CUDACC__
      if ( sigmakin_itry <= maxtry )
        sigmakin_itry++;
      //if ( sigmakin_itry == maxtry )
      //  for (int ihel = 0; ihel < ncomb; ihel++ )
      //    printf( "sigmakin: ihelgood %%2d %%d\n", ihel, sigmakin_goodhel[ihel] );
#endif
      //if ( ighel == 0 ) break; // TEST sectors/requests (issue #16)
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef MGONGPU_CPPSIMD
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      allMEs[ipagV] /= (fptype)denominators;
    }
#else
#ifndef __CUDACC__
    for ( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
      allMEs[ievt] /= (fptype)denominators;
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------
