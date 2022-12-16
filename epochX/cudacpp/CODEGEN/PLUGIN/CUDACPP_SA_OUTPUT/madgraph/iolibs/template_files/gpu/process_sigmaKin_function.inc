    // PART 0 - INITIALISATION (before calculate_wavefunctions)
    // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    allMEs[ievt] = 0;
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    allNumerators[ievt] = 0;
    allDenominators[ievt] = 0;
#endif
#else
    const int npagV = nevt / neppV;
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv = fptype_sv{ 0 };
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
      fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
      fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
      fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
      numerators_sv = fptype_sv{ 0 };
      denominators_sv = fptype_sv{ 0 };
#endif
    }
#endif

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here#ifdef __CUDACC__
#ifdef __CUDACC__
    // *** PART 1a - CUDA (one event per CPU thread) ***
    // Running sum of partial amplitudes squared for event by event color selection (#402)
    // (for the single event or the single SIMD vector of events processed in calculate_wavefunctions)
    fptype_sv jamp2_sv[nParity * ncolor] = { 0 }; // jamp2[nParity][ncolor][neppV] for color choice
    fptype MEs_ighel[ncomb] = { 0 }; // sum of MEs for all good helicities up to ighel (for this event)
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, channelId, allNumerators, allDenominators, jamp2_sv );
#else
      calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, jamp2_sv );
#endif
      MEs_ighel[ighel] = allMEs[ievt];
    }
    // Event-by-event random choice of helicity #403
    //printf( "sigmaKin: ievt=%%4d rndhel=%%f\n", ievt, allrndhel[ievt] );
    for( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      if( allrndhel[ievt] < ( MEs_ighel[ighel] / MEs_ighel[cNGoodHel - 1] ) )
      {
        const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
        allselhel[ievt] = ihelF;
        //printf( "sigmaKin: ievt=%%4d ihel=%%4d\n", ievt, ihelF );
        break;
      }
    }
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    const int channelIdC = channelId - 1; // coloramps.h uses the C array indexing starting at 0
    // Event-by-event random choice of color #402
    fptype targetamp[ncolor] = { 0 };
    for( int icolC = 0; icolC < ncolor; icolC++ )
    {
      if( mgOnGpu::icolamp[channelIdC][icolC] )
        targetamp[icolC] = jamp2_sv[icolC] + ( icolC == 0 ? 0 : targetamp[icolC - 1] );
    }
    //printf( "sigmaKin: ievt=%%4d rndcol=%%f\n", ievt, allrndcol[ievt] );
    for( int icolC = 0; icolC < ncolor; icolC++ )
    {
      if( allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] ) )
      {
        allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
        break;
      }
    }
#endif
#else
    // *** PART 1b - C++ (loop on event pages)
    fptype_sv MEs_ighel[ncomb] = { 0 }; // sum of MEs for all good helicities up to ighel (for the first - and/or only - neppV page)
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Mixed fptypes #537: float for color algebra and double elsewhere
    // Delay color algebra and ME updates (only on even pages)
    assert( npagV %% 2 == 0 );            // SANITY CHECK for mixed fptypes: two neppV-pages are merged to one 2*neppV-page
    const int npagV2 = npagV / 2;        // loop on two SIMD pages (neppV events) at a time
    fptype_sv MEs_ighel2[ncomb] = { 0 }; // sum of MEs for all good helicities up to ighel (for the second neppV page)
#else
    const int npagV2 = npagV;            // loop on one SIMD page (neppV events) at a time
#endif
    /*
#ifdef _OPENMP
    // (NB gcc9 or higher, or clang, is required)
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#pragma omp parallel for default( none ) shared( allmomenta, allcouplings, allMEs, channelId, allNumerators, allDenominators )
#endif // _OPENMP
    */
    for( int ipagV2 = 0; ipagV2 < npagV2; ++ipagV2 )
    {
      // Running sum of partial amplitudes squared for event by event color selection (#402)
      // (for the single event or the single SIMD vector of events processed in calculate_wavefunctions)
      fptype_sv jamp2_sv[nParity * ncolor] = { 0 }; // jamp2[nParity][ncolor][neppV] for color choice
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      const int ievt00 = ipagV2 * neppV * 2; // loop on two SIMD pages (neppV events) at a time
#else
      const int ievt00 = ipagV2 * neppV; // loop on one SIMD page (neppV events) at a time
#endif
      for( int ighel = 0; ighel < cNGoodHel; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
        calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, channelId, allNumerators, allDenominators, jamp2_sv, ievt00 );
#else
        calculate_wavefunctions( ihel, allmomenta, allcouplings, allMEs, jamp2_sv, ievt00 );
#endif
        MEs_ighel[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 ) );
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        MEs_ighel2[ighel] = E_ACCESS::kernelAccess( E_ACCESS::ieventAccessRecord( allMEs, ievt00 + neppV ) );
#endif
      }
      // Event-by-event random choice of helicity #403
      for( int ieppV = 0; ieppV < neppV; ++ieppV )
      {
        const int ievt = ievt00 + ieppV;
        //printf( "sigmaKin: ievt=%%4d rndhel=%%f\n", ievt, allrndhel[ievt] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
#if defined MGONGPU_CPPSIMD
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel][ieppV] / MEs_ighel[cNGoodHel - 1][ieppV] );
#else
          const bool okhel = allrndhel[ievt] < ( MEs_ighel[ighel] / MEs_ighel[cNGoodHel - 1] );
#endif
          if( okhel )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt] = ihelF;
            //printf( "sigmaKin: ievt=%%4d ihel=%%4d\n", ievt, ihelF );
            break;
          }
        }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        const int ievt2 = ievt00 + ieppV + neppV;
        //printf( "sigmaKin: ievt=%%4d rndhel=%%f\n", ievt2, allrndhel[ievt2] );
        for( int ighel = 0; ighel < cNGoodHel; ighel++ )
        {
          if( allrndhel[ievt2] < ( MEs_ighel2[ighel][ieppV] / MEs_ighel2[cNGoodHel - 1][ieppV] ) )
          {
            const int ihelF = cGoodHel[ighel] + 1; // NB Fortran [1,ncomb], cudacpp [0,ncomb-1]
            allselhel[ievt2] = ihelF;
            //printf( "sigmaKin: ievt=%%4d ihel=%%4d\n", ievt, ihelF );
            break;
          }
        }
#endif
      }
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      const int channelIdC = channelId - 1; // coloramps.h uses the C array indexing starting at 0
      // Event-by-event random choice of color #402
      fptype_sv targetamp[ncolor] = { 0 };
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( mgOnGpu::icolamp[channelIdC][icolC] )
          targetamp[icolC] = jamp2_sv[icolC] + ( icolC == 0 ? fptype_sv{ 0 } : targetamp[icolC - 1] );
      }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype_sv targetamp2[ncolor] = { 0 };
      for( int icolC = 0; icolC < ncolor; icolC++ )
      {
        if( mgOnGpu::icolamp[channelIdC][icolC] )
          targetamp2[icolC] = jamp2_sv[ncolor + icolC] + ( icolC == 0 ? fptype_sv{ 0 } : targetamp2[icolC - 1] );
      }
#endif
      for( int ieppV = 0; ieppV < neppV; ++ieppV )
      {
        const int ievt = ievt00 + ieppV;
        //printf( "sigmaKin: ievt=%%4d rndcol=%%f\n", ievt, allrndcol[ievt] );
        for( int icolC = 0; icolC < ncolor; icolC++ )
        {
#if defined MGONGPU_CPPSIMD
          const bool okcol = allrndcol[ievt] < ( targetamp[icolC][ieppV] / targetamp[ncolor - 1][ieppV] );
#else
          const bool okcol = allrndcol[ievt] < ( targetamp[icolC] / targetamp[ncolor - 1] );
#endif
          if( okcol )
          {
            allselcol[ievt] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
            break;
          }
        }
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        const int ievt2 = ievt00 + ieppV + neppV;
        //printf( "sigmaKin: ievt=%%4d rndcol=%%f\n", ievt2, allrndcol[ievt2] );
        for( int icolC = 0; icolC < ncolor; icolC++ )
        {
          if( allrndcol[ievt2] < ( targetamp2[icolC][ieppV] / targetamp2[ncolor - 1][ieppV] ) )
          {
            allselcol[ievt2] = icolC + 1; // NB Fortran [1,ncolor], cudacpp [0,ncolor-1]
            break;
          }
        }
#endif
#endif
      }
    }
#endif

    // PART 2 - FINALISATION (after calculate_wavefunctions)
    // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
#ifdef __CUDACC__
    allMEs[ievt] /= helcolDenominators[0]; // FIXME (#343): assume nprocesses == 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId > 0 ) allMEs[ievt] *= allNumerators[ievt] / allDenominators[ievt]; // FIXME (#343): assume nprocesses == 1
#endif
#else
    for( int ipagV = 0; ipagV < npagV; ++ipagV )
    {
      const int ievt0 = ipagV * neppV;
      fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv /= helcolDenominators[0]; // FIXME (#343): assume nprocesses == 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
      if( channelId > 0 )
      {
        fptype* numerators = NUM_ACCESS::ieventAccessRecord( allNumerators, ievt0 );
        fptype* denominators = DEN_ACCESS::ieventAccessRecord( allDenominators, ievt0 );
        fptype_sv& numerators_sv = NUM_ACCESS::kernelAccess( numerators );
        fptype_sv& denominators_sv = DEN_ACCESS::kernelAccess( denominators );
        MEs_sv *= numerators_sv / denominators_sv; // FIXME (#343): assume nprocesses == 1
      }
#endif
      //for( int ieppV = 0; ieppV < neppV; ieppV++ )
      //{
      //  const unsigned int ievt = ipagV * neppV + ieppV;
      //  printf( "sigmaKin: ievt=%%2d me=%%f\n", ievt, allMEs[ievt] );
      //}
    }
#endif
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------
