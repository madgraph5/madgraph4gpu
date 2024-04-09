! Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
! Created by: J. Alwall (Jul 2010) for the MG5aMC CPP backend.
!==========================================================================
! Copyright (C) 2020-2023 CERN and UCLouvain.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
! Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.
!==========================================================================
! Copyright (C) 2021-2023 Argonne National Laboratory.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: N. Nichols (2021-2023) for the MG5aMC SYCL plugin.
!==========================================================================
      // *** COLOR CHOICE BELOW ***
      // Store the leading color flows for choice of color
      if (jamp2_sv) { // disable color choice if nullptr
          for (size_t icolC = 0; icolC < ncolor; icolC++) {
              jamp2_sv[icolC] += CXABS2(jamp_sv[icolC]);
          }
      }
      // *** COLOR MATRIX BELOW ***
      // (This method used to be called %(process_class_name)s::matrix_%(proc_name)s(%(matrix_args)s)?)

      // The color matrix (initialize all array elements, with ncolor=1)
      // [NB do keep 'static' for these constexpr arrays, see issue #283]
%(color_matrix_lines)s

      // Sum and square the color flows to get the matrix element
      // (compute |M|^2 by squaring |M|, taking into account colours)
      fptype_sv deltaMEs = FPZERO_SV;
      for( size_t icol = 0; icol < ncolor; icol++ ) {
          cxtype_sv ztemp_sv = CXZERO_SV;
          for( size_t jcol = 0; jcol < ncolor; jcol++ ) {
              ztemp_sv += cf[icol][jcol]*jamp_sv[jcol];
          }
          // OLD implementation: why is this not slower? maybe compiler does not compute imaginary part of "ztemp_sv*cxconj(jamp_sv[icol])"?
          //deltaMEs += cxreal( ztemp_sv * cxconj( jamp_sv[icol] ) ) / denom[icol];
          // NEW implementation: keep this even if (surprisingly) it is not faster! it is clearer and easier for tensor core offload anyway...
          // Rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB!
          deltaMEs += (CXREAL(ztemp_sv)*CXREAL(jamp_sv[icol]) + CXIMAG(ztemp_sv)*CXIMAG(jamp_sv[icol]))/denom[icol];
      }
      //FIXME test if faster
      //for (size_t icol = 0; icol < ncolor; icol++) {
      //    fptype_sv ztempR_sv = FPZERO_SV;
      //    fptype_sv ztempI_sv = FPZERO_SV;
      //    for (size_t jcol = 0; jcol < ncolor; jcol++ ) {
      //        fptype_sv jampRj_sv = CXREAL(jamp_sv[jcol]);
      //        fptype_sv jampIj_sv = CXIMAG(jamp_sv[jcol]);
      //        ztempR_sv += cf[icol][jcol]*jampRj_sv;
      //        ztempI_sv += cf[icol][jcol]*jampIj_sv;
      //    }
      //    deltaMEs += ( ztempR_sv*CXREAL(jamp_sv[icol]) + ztempI_sv*CXIMAG(jamp_sv[icol]))/denom[icol];
      //}

      // *** STORE THE RESULTS ***

      // Store the leading color flows for choice of color
      // (NB: jamp2_sv must be an array of fptype_sv)
      // for( int icol = 0; icol < ncolor; icol++ )
      // jamp2_sv[0][icol] += cxreal( jamp_sv[icol]*cxconj( jamp_sv[icol] ) );

      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event(s)
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs += deltaMEs;
      return allMEs;
  }
