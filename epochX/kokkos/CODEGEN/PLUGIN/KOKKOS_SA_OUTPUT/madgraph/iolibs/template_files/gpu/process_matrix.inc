
  // *** COLOR ALGEBRA BELOW ***
  // (This method used to be called %(process_class_name)s::matrix_%(proc_name)s(%(matrix_args)s)?)

  // The color matrix
%(color_matrix_lines)s

  // Sum and square the color flows to get the matrix element
  // (compute |M|^2 by squaring |M|, taking into account colours)
  fptype deltaMEs = { 0 }; // all zeros
  for( int icol = 0; icol < ncolor; icol++ )
  {
    cxtype ztemp = cxmake( 0, 0 );
    for( int jcol = 0; jcol < ncolor; jcol++ )
      ztemp += cf[icol][jcol] * jamp[jcol];
    deltaMEs += cxreal( ztemp * cxconj( jamp[icol] ) ) / denom[icol];
  }

  // *** STORE THE RESULTS ***

  // Store the leading color flows for choice of color
  // (NB: jamp2 must be an array of fptype)
  // for( int icol = 0; icol < ncolor; icol++ )
  // jamp2[0][icol] += cxreal( jamp[icol]*cxconj( jamp[icol] ) );

  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
  // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
  allMEs += deltaMEs;
  //printf( "calculate_wavefunction: %%6d %%2d %%f\n", ievt, ihel, allMEs[ievt] );


  return allMEs;
}
