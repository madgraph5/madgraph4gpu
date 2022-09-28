    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    const int nevt = gpublocks * gputhreads;
    cudaMemset( allMEs, 0, nevt * sizeof( fptype ) );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    cudaMemset( allNumerators, 0, nevt * sizeof( fptype ) );
    cudaMemset( allDenominators, 0, nevt * sizeof( fptype ) );
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
      calculate_wavefunctions<<<gpublocks, gputhreads>>>( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, channelId );
#else
      calculate_wavefunctions<<<gpublocks, gputhreads>>>( ihel, allmomenta, allcouplings, allMEs );
#endif
#else
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, allNumerators, allDenominators, channelId, nevt );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, nevt );
#endif
#endif
      //if ( ighel == 0 ) break; // TEST sectors/requests (issue #16)
    }

    // TODO-OM: See how to handle the renormalization here...
    // Dedicated kernel or move it within each calculate_wavefunctions?

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    normalise_output<<<gpublocks, gputhreads>>>( allMEs, allNumerators, allDenominators, channelId, denominators[0] );
#else
    normalise_output<<<gpublocks, gputhreads>>>( allMEs, denominators[0] );
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
#endif
        //printf( "sigmaKin: ievt=%%2d me=%%f\n", ievt, allMEs[ievt] );
      }
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__ /* clang-format off */
  __global__ void
  normalise_output( fptype* allMEs,                // output: allMEs[nevt], |M|^2 running_sum_over_helicities
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
                    const fptype* allNumerators,   // output: multichannel numerators[nevt], running_sum_over_helicities
                    const fptype* allDenominators, // output: multichannel denominators[nevt], running_sum_over_helicities
                    const unsigned int channelId,  // input: multichannel channel id (1 to #diagrams); 0 to disable channel enhancement
#endif
                    const fptype globaldenom ) /* clang-format on */
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread)
    allMEs[ievt] /= globaldenom;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // FIXME (#343): assume nprocesses == 1
#endif
    return;
  }
#endif

  //--------------------------------------------------------------------------
