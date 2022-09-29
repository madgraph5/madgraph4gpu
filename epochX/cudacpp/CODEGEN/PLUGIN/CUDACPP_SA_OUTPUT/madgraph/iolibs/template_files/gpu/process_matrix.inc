
      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called %(process_class_name)s::matrix_%(proc_name)s(%(matrix_args)s)?)
%(color_matrix_lines)s

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = { 0 }; // all zeros https://en.cppreference.com/w/c/language/array_initialization#Notes
      for( int icol = 0; icol < ncolor; icol++ )
      {
        fptype_sv ztempR_sv = cf[icol][icol] * jampR_sv[icol];
        fptype_sv ztempI_sv = cf[icol][icol] * jampR_sv[icol];
	
        for( int jcol = icol+1; jcol < ncolor; jcol++ ){
	  // factor of 2 here because we only loop over the up diagonal part of the matrix
          ztempR_sv += 2*cf[icol][jcol] * jampR_sv[jcol];
	  ztempI_sv += 2*cf[icol][jcol] * jampI_sv[jcol];
        }
	// here we use the property that M is a real matrix such that we can
        // Rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB!
	deltaMEs += (jampR_sv[icol]*ztempR_sv + jampI_sv[icol]*ztempI_sv)/denom[icol];
	}
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
