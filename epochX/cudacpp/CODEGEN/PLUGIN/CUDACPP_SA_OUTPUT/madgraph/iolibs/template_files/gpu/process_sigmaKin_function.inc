#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    for( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "sigmakin: ievt %%d\n", ievt );
#endif

      // Denominators: spins, colors and identical particles
      const int denominators = %(den_factors)s; // FIXME: assume process.nprocesses == 1 for the moment (eventually denominators[nprocesses]?)

      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype allMEsThisEvt = 0;

#ifdef __CUDACC__
      // CUDA - using precomputed good helicities
      for ( int ighel = 0; ighel < cNGoodHel; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        calculate_wavefunctions( ihel, allmomenta, allMEsThisEvt );
      }
#else
      // C++ - compute good helicities within this loop
      fptype allMEsThisEvtLast = 0; // check for good helicities
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry > maxtry && !sigmakin_goodhel[ihel] ) continue;
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
        calculate_wavefunctions( ihel, allmomenta, allMEsThisEvt, ievt );
        if ( sigmakin_itry <= maxtry )
        {
          if ( !sigmakin_goodhel[ihel] && allMEsThisEvt > allMEsThisEvtLast ) sigmakin_goodhel[ihel] = true;
          allMEsThisEvtLast = allMEsThisEvt;
        }
      }
#endif

      // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEsThisEvt /= denominators;

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
      // ** END LOOP ON IEVT **
    }
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------
