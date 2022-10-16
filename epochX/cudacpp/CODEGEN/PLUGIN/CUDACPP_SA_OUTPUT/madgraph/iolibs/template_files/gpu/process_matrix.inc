
      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called %(process_class_name)s::matrix_%(proc_name)s(%(matrix_args)s)?)
%(color_matrix_lines)s

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
      // Use the property that M is a real matrix (see #475):
      // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
      for( int icol = 0; icol < ncolor; icol++ )
      {
        // In addition, for C++ use the property that M is symmetric (see #475):
        // we gain (a factor 2?) in speed here as we only loop over the up diagonal part of the matrix!
        // For CUDA, the old implementation is (surprisingly?!) faster expecially for complex physics and large grids
#ifndef __CUDACC__
        // Diagonal terms
        fptype_sv jampRi_sv = cxreal( jamp_sv[icol] );
        fptype_sv jampIi_sv = cximag( jamp_sv[icol] );
        fptype_sv ztempR_sv = cf[icol][icol] * jampRi_sv;
        fptype_sv ztempI_sv = cf[icol][icol] * jampIi_sv;
        // Off-diagonal terms
        for( int jcol = icol + 1; jcol < ncolor; jcol++ )
        {
          fptype_sv jampRj_sv = cxreal( jamp_sv[jcol] );
          fptype_sv jampIj_sv = cximag( jamp_sv[jcol] );
          ztempR_sv += 2 * cf[icol][jcol] * jampRj_sv;
          ztempI_sv += 2 * cf[icol][jcol] * jampIj_sv;
        }
        deltaMEs += ( jampRi_sv * ztempR_sv + jampIi_sv * ztempI_sv ) / denom[icol];
#else
        cxtype_sv ztemp_sv = cxzero_sv();
        for( int jcol = 0; jcol < ncolor; jcol++ )
          ztemp_sv += cf[icol][jcol] * jamp_sv[jcol];
        deltaMEs += ( cxreal( ztemp_sv ) * cxreal( jamp_sv[icol] ) + cximag( ztemp_sv ) * cximag( jamp_sv[icol] ) ) / denom[icol];
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
