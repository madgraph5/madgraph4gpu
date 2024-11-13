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

  // *** COLOR MATRIX BELOW ***
%(color_matrix_lines)s

#ifdef MGONGPUCPP_GPUIMPL
  // The normalized color matrix (divide each column by denom)
  template<typename T>
  struct NormalizedColorMatrix
  {
    constexpr __device__ NormalizedColorMatrix()
      : value()
    {
      for( int icol = 0; icol < ncolor; icol++ )
        for( int jcol = 0; jcol < ncolor; jcol++ )
          value[icol * ncolor + jcol] = colorMatrix[icol][jcol] / colorDenom[icol];
    }
    T value[ncolor * ncolor];
  };
#ifndef MGONGPU_HAS_NO_BLAS
  // The fptype version is only used by BLAS (TEMPORARY: mixed mode is not fully supported and uses double precision for BLAS color sums)
  static __device__ fptype s_pNormalizedColorMatrix[ncolor * ncolor];
#endif
  // The fptype2 version is the default used by kernels (supporting mixed floating point mode)
  static __device__ fptype2 s_pNormalizedColorMatrix2[ncolor * ncolor];
#endif

  //--------------------------------------------------------------------------

#ifndef MGONGPUCPP_GPUIMPL
  INLINE void
  color_sum_cpu( fptype* allMEs,              // output: allMEs[nevt], add |M|^2 for this specific helicity
                 const cxtype_sv* allJamp_sv, // input: jamp_sv[ncolor] (float/double) or jamp_sv[2*ncolor] (mixed) for one specific helicity
                 const int ievt0 )            // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
  {
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
          value[icol][icol] = colorMatrix[icol][icol] / colorDenom[icol];
          // Off-diagonal terms
          for( int jcol = icol + 1; jcol < ncolor; jcol++ )
            value[icol][jcol] = 2 * colorMatrix[icol][jcol] / colorDenom[icol];
        }
      }
      fptype2 value[ncolor][ncolor];
    };
    static constexpr auto cf2 = TriangularNormalizedColorMatrix();
    // Use the property that M is a real matrix (see #475):
    // we can rewrite the quadratic form (A-iB)(M)(A+iB) as AMA - iBMA + iBMA + BMB = AMA + BMB
    // In addition, on C++ use the property that M is symmetric (see #475),
    // and also use constexpr to compute "2*" and "/colorDenom[icol]" once and for all at compile time:
    // we gain (not a factor 2...) in speed here as we only loop over the up diagonal part of the matrix.
    // Strangely, CUDA is slower instead, so keep the old implementation for the moment.
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
    // *** STORE THE RESULTS ***
    using E_ACCESS = HostAccessMatrixElements; // non-trivial access: buffer includes all events
    fptype* MEs = E_ACCESS::ieventAccessRecord( allMEs, ievt0 );
    // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
    fptype_sv& MEs_sv = E_ACCESS::kernelAccess( MEs );
    MEs_sv += deltaMEs; // fix #435
#if defined MGONGPU_CPPSIMD and defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    fptype* MEs_next = E_ACCESS::ieventAccessRecord( allMEs, ievt0 + neppV );
    fptype_sv& MEs_sv_next = E_ACCESS::kernelAccess( MEs_next );
    MEs_sv_next += deltaMEs_next;
#endif
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  __global__ INLINE void
  color_sum_kernel( fptype* allMEs,          // output: allMEs[nevt], add |M|^2 for this specific helicity
                    const fptype* allJamps ) // input: jamp[ncolor*2*nevt] for one specific helicity
  {
    using J_ACCESS = DeviceAccessJamp;
    fptype jampR[ncolor];
    fptype jampI[ncolor];
    for( int icol = 0; icol < ncolor; icol++ )
    {
      cxtype jamp = J_ACCESS::kernelAccessIcolConst( allJamps, icol );
      jampR[icol] = jamp.real();
      jampI[icol] = jamp.imag();
    }
    // Loop over icol
    fptype deltaMEs = { 0 };
    for( int icol = 0; icol < ncolor; icol++ )
    {
      fptype2 ztempR = { 0 };
      fptype2 ztempI = { 0 };
      // Loop over jcol
      for( int jcol = 0; jcol < ncolor; jcol++ )
      {
        fptype2 jampRj = jampR[jcol];
        fptype2 jampIj = jampI[jcol];
        ztempR += s_pNormalizedColorMatrix2[icol * ncolor + jcol] * jampRj; // use fptype2 version of color matrix
        ztempI += s_pNormalizedColorMatrix2[icol * ncolor + jcol] * jampIj; // use fptype2 version of color matrix
      }
      deltaMEs += ztempR * jampR[icol];
      deltaMEs += ztempI * jampI[icol];
    }
    // *** STORE THE RESULTS ***
    using E_ACCESS = DeviceAccessMatrixElements; // non-trivial access: buffer includes all events
    // NB: color_sum ADDS |M|^2 for one helicity to the running sum of |M|^2 over helicities for the given event(s)
    E_ACCESS::kernelAccess( allMEs ) += deltaMEs; // fix #435
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
#ifndef MGONGPU_HAS_NO_BLAS
  INLINE void
  color_sum_blas( fptype* allMEs,               // output: allMEs[nevt], add |M|^2 for this specific helicity
                  const fptype* allJamps,       // input: jamp[ncolor*2*nevt] for one specific helicity
                  fptype* allBlasTmp,           // tmp: blasTmp[ncolor*2*nevt] for one specific helicity
                  gpuBlasHandle_t* pBlasHandle, // input: cuBLAS/hipBLAS handle
                  const int gpublocks,          // input: cuda gpublocks
                  const int gputhreads )        // input: cuda gputhreads
  {
    const int nevt = gpublocks * gputhreads;

    // Get the address associated with the normalized color matrix in device memory
    static fptype* devNormColMat = nullptr;
    if( !devNormColMat ) gpuGetSymbolAddress( (void**)&devNormColMat, s_pNormalizedColorMatrix );

    // New striding for cuBLAS from DeviceAccessJamp:
    // - allJamps(icol,ievt).real is allJamps[0 * ncolor * nevt + icol * nevt + ievt]
    // - allJamps(icol,ievt).imag is allJamps[1 * ncolor * nevt + icol * nevt + ievt]
    const fptype* allJampsReal = allJamps;
    const fptype* allJampsImag = allJamps + ncolor * nevt;
    fptype* allTmpReal = allBlasTmp;
    fptype* allTmpImag = allBlasTmp + ncolor * nevt;

    // Step 1: Compute Tmp[ncolor][nevt] = ColorMatrix[ncolor][ncolor] * JampsVector[ncolor][nevt] for both real and imag
    // In this case alpha=1 and beta=0 (the operation is Tmp = alpha * ColorMatrix * JampsVector + beta * Tmp)
    fptype alpha1 = 1;
    fptype beta1 = 0;
    const int ncolorM = ncolor;
    const int nevtN = nevt;
    const int ncolorK = ncolor;
    checkGpuBlas( gpuBlasTgemm( *pBlasHandle,
                                CUBLAS_OP_N,             // do not transpose ColMat
                                CUBLAS_OP_N,             // do not transpose JampsV
                                ncolorM, nevtN, ncolorK,
                                &alpha1,
                                devNormColMat, ncolorM,  // ColMat is ncolorM x ncolorK
                                allJampsReal, ncolorK,   // JampsV is ncolorK x nevtN
                                &beta1,
                                allTmpReal, ncolorM ) ); // Tmp is ncolorM x ncolorN
    checkGpuBlas( gpuBlasTgemm( *pBlasHandle,
                                CUBLAS_OP_N,             // do not transpose ColMat
                                CUBLAS_OP_N,             // do not transpose JampsV
                                ncolorM, nevtN, ncolorK,
                                &alpha1,
                                devNormColMat, ncolorM,  // ColMat is ncolorM x ncolorK
                                allJampsImag, ncolorK,   // JampsV is ncolorK x nevtN
                                &beta1,
                                allTmpImag, ncolorM ) ); // Tmp is ncolorM x ncolorN

    // Step 2: For each ievt, compute the dot product of JampsVector[ncolor][ievt] dot tmp[ncolor][ievt]
    // In this case alpha=1 and beta=1 (the operation is ME = alpha * ( Tmp dot JampsVector ) + beta * ME)
    // Use cublasSgemmStridedBatched to perform these batched dot products in one call
    fptype alpha2 = 1;
    fptype beta2 = 1;
    checkGpuBlas( gpuBlasTgemmStridedBatched( *pBlasHandle,
                                              CUBLAS_OP_T,                  // transpose JampsV
                                              CUBLAS_OP_N,                  // do not transpose Tmp
                                              1, 1, ncolor,                 // result is 1x1 (dot product)
                                              &alpha2,
                                              allJampsReal, ncolor, ncolor, // allJamps is ncolor x nevt, with stride ncolor for each ievt column
                                              allTmpReal, ncolor, ncolor,   // allTmp is ncolor x nevt, with stride ncolor for each ievt column
                                              &beta2,
                                              allMEs, 1, 1,                 // output is a 1x1 result for each "batch" (i.e. for each ievt)
                                              nevt ) );                     // there are nevt "batches"
    checkGpuBlas( gpuBlasTgemmStridedBatched( *pBlasHandle,
                                              CUBLAS_OP_T,                  // transpose JampsV
                                              CUBLAS_OP_N,                  // do not transpose Tmp
                                              1, 1, ncolor,                 // result is 1x1 (dot product)
                                              &alpha2,
                                              allJampsImag, ncolor, ncolor, // allJamps is ncolor x nevt, with stride ncolor for each ievt column
                                              allTmpImag, ncolor, ncolor,   // allTmp is ncolor x nevt, with stride ncolor for each ievt column
                                              &beta2,
                                              allMEs, 1, 1,                 // output is a 1x1 result for each "batch" (i.e. for each ievt)
                                              nevt ) );                     // there are nevt "batches"

  }
#endif /* clang-format on */
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  INLINE void
  color_sum_gpu( fptype* allMEs,               // output: allMEs[nevt], add |M|^2 for this specific helicity
                 const fptype* allJamps,       // input: jamp[ncolor*2*nevt] for one specific helicity
                 fptype* allBlasTmp,           // tmp: blasTmp[ncolor*2*nevt] for one specific helicity
                 gpuStream_t stream,           // input: cuda stream (nullptr indicates the default stream)
                 gpuBlasHandle_t* pBlasHandle, // input: cuBLAS/hipBLAS handle
                 const int gpublocks,          // input: cuda gpublocks
                 const int gputhreads )        // input: cuda gputhreads
  {
#ifdef MGONGPU_HAS_NO_BLAS
    assert( allBlasTmp == nullptr );  // sanity check for HASBLAS=hasNoBlas
    assert( pBlasHandle == nullptr ); // sanity check for HASBLAS=hasNoBlas
#endif
    if( !pBlasHandle ) // HASBLAS=hasNoBlas or CUDACPP_RUNTIME_BLASCOLORSUM not set
    {
      assert( allBlasTmp == nullptr );
      gpuLaunchKernelStream( color_sum_kernel, gpublocks, gputhreads, stream, allMEs, allJamps );
    }
#ifndef MGONGPU_HAS_NO_BLAS
    else
    {
      assert( allBlasTmp != nullptr );
      color_sum_blas( allMEs, allJamps, allBlasTmp, pBlasHandle, gpublocks, gputhreads );
    }
#endif
  }
#endif
