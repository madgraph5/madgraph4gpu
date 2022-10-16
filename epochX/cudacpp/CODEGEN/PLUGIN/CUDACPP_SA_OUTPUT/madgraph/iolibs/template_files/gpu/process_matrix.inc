
      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called %(process_class_name)s::matrix_%(proc_name)s(%(matrix_args)s)?)
%(color_matrix_lines)s

#ifndef __CUDACC__
      // Pre-compute a constexpr triangular color matrix properly normalized #475
      struct TriangularNormalizedColorMatrix
      {
        // See https://stackoverflow.com/a/34465458
        __host__ __device__ constexpr TriangularNormalizedColorMatrix()
          : value()
        {
          for( int icol = 0; icol < ncolor; icol++ )
          {
            // Diagonal terms
            value[icol][icol] = cf[icol][icol] / denom[icol];
            // Off-diagonal terms
            for( int jcol = icol + 1; jcol < ncolor; jcol++ )
              value[icol][jcol] = 2 * cf[icol][jcol] / denom[icol];
          }
        }
        fptype value[ncolor][ncolor];
      };
      static constexpr auto cf2 = TriangularNormalizedColorMatrix();
#endif

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
      // Use the property that M is a real matrix (see #475):
      // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
      // In addition, on C++ use the property that M is symmetric (see #475),
      // and also use constexpr to compute "2*" and "/denom[icol]" once and for all at compile time:
      // we gain (not a factor 2...) in speed here as we only loop over the up diagonal part of the matrix.
      // Strangely, CUDA is slower instead, so keep the old implementation for the moment.
      for( int icol = 0; icol < ncolor; icol++ )
      {
#ifndef __CUDACC__
        // Diagonal terms
        fptype_sv jampRi_sv = cxreal( jamp_sv[icol] );
        fptype_sv jampIi_sv = cximag( jamp_sv[icol] );
        fptype_sv ztempR_sv = cf2.value[icol][icol] * jampRi_sv;
        fptype_sv ztempI_sv = cf2.value[icol][icol] * jampIi_sv;
        // Off-diagonal terms
        for( int jcol = icol + 1; jcol < ncolor; jcol++ )
        {
          fptype_sv jampRj_sv = cxreal( jamp_sv[jcol] );
          fptype_sv jampIj_sv = cximag( jamp_sv[jcol] );
          ztempR_sv += cf2.value[icol][jcol] * jampRj_sv;
          ztempI_sv += cf2.value[icol][jcol] * jampIj_sv;
        }
        deltaMEs += ( ztempR_sv * jampRi_sv + ztempI_sv * jampIi_sv );
#else
        fptype_sv ztempR_sv = { 0 };
        fptype_sv ztempI_sv = { 0 };
        for( int jcol = 0; jcol < ncolor; jcol++ )
        {
          ztempR_sv += cf[icol][jcol] * cxreal( jamp_sv[jcol] );
          ztempI_sv += cf[icol][jcol] * cximag( jamp_sv[jcol] );
        }
        deltaMEs += ( ztempR_sv * cxreal( jamp_sv[icol] ) + ztempI_sv * cximag( jamp_sv[icol] ) ) / denom[icol];
#endif
      }

      // *** STORE THE RESULTS ***

      // Store the leading color flows for choice of color
      // (NB: jamp2_sv must be an array of fptype_sv)
      // for( int icol = 0; icol < ncolor; icol++ )
      // jamp2_sv[0][icol] += cxreal( jamp_sv[icol]*cxconj( jamp_sv[icol] ) );

      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
      MEs_sv += deltaMEs; // fix #435
      /*
#ifdef __CUDACC__
      if ( cNGoodHel > 0 ) printf( "calculate_wavefunctions: ievt=%%6d ihel=%%2d me_running=%%f\n", blockDim.x * blockIdx.x + threadIdx.x, ihel, MEs_sv );
#else
#ifdef MGONGPU_CPPSIMD
      if( cNGoodHel > 0 )
        for( int ieppV = 0; ieppV < neppV; ieppV++ )
          printf( "calculate_wavefunctions: ievt=%%6d ihel=%%2d me_running=%%f\n", ipagV * neppV + ieppV, ihel, MEs_sv[ieppV] );
#else
      if ( cNGoodHel > 0 ) printf( "calculate_wavefunctions: ievt=%%6d ihel=%%2d me_running=%%f\n", ipagV, ihel, MEs_sv );
#endif
#endif
      */
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }
