    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    allMEs[ievt] = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    allNumerators[ievt] = 0;
    allDenominators[ievt] = 0;
    for( int icol = 0; icol < ncolor; ++icol )
    {
	jamp2[ievt][icol] = 0;
	}		    
#endif
#else
    const int npagV = nevt / neppV;
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      for( int ieppV = 0; ieppV < neppV; ieppV++ )
      {
        const unsigned int ievt = ipagV * neppV + ieppV;
        allMEs[ievt] = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        allNumerators[ievt] = 0;
        allDenominators[ievt] = 0;
	for( int icol = 0; icol < ncolor; ++icol )
    	{	
		jamp2[ievt][icol] = 0;
	}		    
#endif
      }
    }
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef __CUDACC__
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, jamp2, channelId );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs );
#endif
#else
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, jamp2, channelId, nevt );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, nevt );
#endif
#endif
      //if ( ighel == 0 ) break; // TEST sectors/requests (issue #16)
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    allMEs[ievt] /= denominators[0]; // FIXME (#343): assume nprocesses == 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // FIXME (#343): assume nprocesses == 1
#endif
#else
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      for( int ieppV = 0; ieppV < neppV; ieppV++ )
      {
        const unsigned int ievt = ipagV * neppV + ieppV;
        allMEs[ievt] /= denominators[0];                                                 // FIXME (#343): assume nprocesses == 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        if( channelId > 0 ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // FIXME (#343): assume nprocesses == 1
	fptype is_LC[%(nb_channel)s][%(nb_color)s] = %(is_LC)s ;

	int nb_color = %(nb_color)i;
	float sumjamp2 = 0.;
	for (int c=0; c < nb_color; c++){
	if (is_LC[channelId][c]){                       // FIXME (#343): assume nprocesses == 1
		sumjamp2 += jamp2[c][ievt];
	}
	}
	sumjamp2 = sumjamp2*random_for_color[ievt];
	for (int c=0; c < nb_color; c++){
		if (is_LC[channelId][c]){
		   sumjamp2 -= jamp2[c][ievt];
		   if (sumjamp2 < 0){
		      selected_color[ievt] = c;
		      break;
		      }
		 }
	}
#endif
        //printf( "sigmaKin: ievt=%%2d me=%%f\n", ievt, allMEs[ievt] );
      }
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------
