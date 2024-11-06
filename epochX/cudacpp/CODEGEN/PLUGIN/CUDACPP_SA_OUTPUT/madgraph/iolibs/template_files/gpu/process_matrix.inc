! Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
! Created by: J. Alwall (Jul 2010) for the MG5aMC CPP backend.
!==========================================================================
! Copyright (C) 2020-2024 CERN and UCLouvain.
! Licensed under the GNU Lesser General Public License (version 3 or later).
! Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
! Further modified by: J. Teig, A. Valassi (2021-2024) for the MG5aMC CUDACPP plugin.
!==========================================================================

      // *** COLOR CHOICE BELOW ***

      // Store the leading color flows for choice of color
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#ifndef MGONGPUCPP_GPUIMPL
      if( jamp2_sv ) // disable color choice if nullptr
      {
        for( int icol = 0; icol < ncolor; icol++ )
          jamp2_sv[ncolor * iParity + icol] += cxabs2( jamp_sv[icol] ); // may underflow #831
      }
#else
      assert( iParity == 0 ); // sanity check for J2_ACCESS
      using J2_ACCESS = DeviceAccessJamp2;
      if( colAllJamp2s ) // disable color choice if nullptr
      {
        for( int icol = 0; icol < ncolor; icol++ )
          // NB: atomicAdd is needed after moving to cuda streams with one helicity per stream!
          atomicAdd( &J2_ACCESS::kernelAccessIcol( colAllJamp2s, icol ), cxabs2( jamp_sv[icol] ) );
      }
#endif
#endif

      // *** PREPARE OUTPUT JAMPS ***
#ifdef MGONGPUCPP_GPUIMPL
      // In CUDA, copy the local jamp to the output global-memory jamp
      using J_ACCESS = DeviceAccessJamp;
      for( int icol = 0; icol < ncolor; icol++ )
        J_ACCESS::kernelAccessIcol( allJamps, icol ) = jamp_sv[icol];
#else
      // In C++, copy the local jamp to the output array passed as function argument
      for( int icol = 0; icol < ncolor; icol++ )
        allJamp_sv[iParity * ncolor + icol] = jamp_sv[icol];
#endif
    }
    // END LOOP ON IPARITY

    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __global__ INLINE void /* clang-format off */
  color_sum( fptype* allMEs,              // output: allMEs[nevt], add |M|^2 for this specific helicity
#ifdef MGONGPUCPP_GPUIMPL
             const fptype* allJamps       // input: jamp[ncolor*2*nevt] for one specific helicity
#else
             const cxtype_sv* allJamp_sv, // input: jamp_sv[ncolor] (float/double) or jamp_sv[2*ncolor] (mixed) for one specific helicity
             const int ievt0              // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
#endif
             ) /* clang-format on */
  {
#ifdef MGONGPUCPP_GPUIMPL
    //using namespace mg5amcGpu;
    using E_ACCESS = DeviceAccessMatrixElements; // non-trivial access: buffer includes all events
#else
    //using namespace mg5amcCpu;
    using E_ACCESS = HostAccessMatrixElements; // non-trivial access: buffer includes all events
#endif

    // *** COLOR MATRIX BELOW ***
%(color_matrix_lines)s

#ifndef MGONGPUCPP_GPUIMPL
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
      fptype2 value[ncolor][ncolor];
    };
    static constexpr auto cf2 = TriangularNormalizedColorMatrix();
#endif

    // Use the property that M is a real matrix (see #475):
    // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
    // In addition, on C++ use the property that M is symmetric (see #475),
    // and also use constexpr to compute "2*" and "/denom[icol]" once and for all at compile time:
    // we gain (not a factor 2...) in speed here as we only loop over the up diagonal part of the matrix.
    // Strangely, CUDA is slower instead, so keep the old implementation for the moment.

#ifndef MGONGPUCPP_GPUIMPL

    // === C++ START ===
    fptype_sv deltaMEs = { 0 };
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    fptype_sv deltaMEs_next = { 0 };
    // Mixed mode: merge two neppV vectors into one neppV2 vector
    fptype2_sv jampR_sv[ncolor];
    fptype2_sv jampI_sv[ncolor];
    for( int icol = 0; icol < ncolor; icol++ )
    {
      jampR_sv[icol] = fpvmerge( cxreal( allJamp_sv[icol] ), cxreal( allJamp_sv[ncolor + icol] ) );
      jampI_sv[icol] = fpvmerge( cximag( allJamp_sv[icol] ), cximag( allJamp_sv[ncolor + icol] ) );
    }
#else
    const cxtype_sv* jamp_sv = allJamp_sv;
#endif
    // Loop over icol
    for( int icol = 0; icol < ncolor; icol++ )
    {
      // Diagonal terms
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      fptype2_sv& jampRi_sv = jampR_sv[icol];
      fptype2_sv& jampIi_sv = jampI_sv[icol];
#else
      fptype2_sv jampRi_sv = (fptype2_sv)( cxreal( jamp_sv[icol] ) );
      fptype2_sv jampIi_sv = (fptype2_sv)( cximag( jamp_sv[icol] ) );
#endif
      fptype2_sv ztempR_sv = cf2.value[icol][icol] * jampRi_sv;
      fptype2_sv ztempI_sv = cf2.value[icol][icol] * jampIi_sv;
      // Loop over jcol
      for( int jcol = icol + 1; jcol < ncolor; jcol++ )
      {
        // Off-diagonal terms
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
        fptype2_sv& jampRj_sv = jampR_sv[jcol];
        fptype2_sv& jampIj_sv = jampI_sv[jcol];
#else
        fptype2_sv jampRj_sv = (fptype2_sv)( cxreal( jamp_sv[jcol] ) );
        fptype2_sv jampIj_sv = (fptype2_sv)( cximag( jamp_sv[jcol] ) );
#endif
        ztempR_sv += cf2.value[icol][jcol] * jampRj_sv;
        ztempI_sv += cf2.value[icol][jcol] * jampIj_sv;
      }
      fptype2_sv deltaMEs2 = ( jampRi_sv * ztempR_sv + jampIi_sv * ztempI_sv ); // may underflow #831
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
      deltaMEs += fpvsplit0( deltaMEs2 );
      deltaMEs_next += fpvsplit1( deltaMEs2 );
#else
      deltaMEs += deltaMEs2;
#endif
    }
    // === C++ END ===

#else

    // === CUDA START ===
    fptype_sv deltaMEs = { 0 };
    using J_ACCESS = DeviceAccessJamp;
    cxtype jamp_sv[ncolor];
    for( int icol = 0; icol < ncolor; icol++ )
      jamp_sv[icol] = J_ACCESS::kernelAccessIcolConst( allJamps, icol );
    // Loop over icol
    for( int icol = 0; icol < ncolor; icol++ )
    {
      fptype2_sv ztempR_sv = { 0 };
      fptype2_sv ztempI_sv = { 0 };
      // Loop over jcol
      for( int jcol = 0; jcol < ncolor; jcol++ )
      {
        fptype2_sv jampRj_sv = cxreal( jamp_sv[jcol] );
        fptype2_sv jampIj_sv = cximag( jamp_sv[jcol] );
        ztempR_sv += cf[icol][jcol] * jampRj_sv;
        ztempI_sv += cf[icol][jcol] * jampIj_sv;
      }
      deltaMEs += ( ztempR_sv * cxreal( jamp_sv[icol] ) + ztempI_sv * cximag( jamp_sv[icol] ) ) / denom[icol];
    }
    // === CUDA END ===

#endif

    // *** STORE THE RESULTS ***
#ifndef MGONGPUCPP_GPUIMPL
    fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
#else
    fptype* MEs = allMEs;
#endif
    // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
    fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
    MEs_sv += deltaMEs; // fix #435
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    fptype* MEs_next = E_ACCESS::ieventAccessRecord( allMEs, ievt0 + neppV );
    fptype_sv& MEs_sv_next = E_ACCESS::kernelAccess( MEs_next );
    MEs_sv_next += deltaMEs_next;
#endif
  }
