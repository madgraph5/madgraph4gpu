
      // *** COLOR ALGEBRA BELOW ***
      // (This method used to be called %(process_class_name)s::matrix_%(proc_name)s(%(matrix_args)s)?)

      // The color matrix
%(color_matrix_lines)s

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      for( int icol=0; icol < ncolor; icol++ )
      {
        cxtype ztemp = cxmake( 0, 0 );
        for( int jcol = 0; jcol < ncolor; jcol++ )
          ztemp = ztemp + cf[icol][jcol] * jamp_sv[jcol];
        allMEs = allMEs + cxreal(ztemp*conj(jamp_sv[icol])) / denom[icol];
      }

      // Store the leading color flows for choice of color
      //for( i=0; i < ncolor; i++ )
      //  jamp2_sv[%(proc_number)d][i] += real(jamp_sv[i]*conj(jamp_sv[i]));

    }

    mgDebug( 1, __FUNCTION__ );
    return;
  }
