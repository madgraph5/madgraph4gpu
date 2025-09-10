// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "color_sum.h"

#include "mgOnGpuConfig.h"

#include "MemoryAccessMatrixElements.h"

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors

  //--------------------------------------------------------------------------

  // *** COLOR MATRIX BELOW ***

  // The color denominators (initialize all array elements, with ncolor=1)
  // [NB do keep 'static' for these constexpr arrays, see issue #283]
  static constexpr fptype2 colorDenom[ncolor] = { 1 }; // 1-D array[1]

  // The color matrix (initialize all array elements, with ncolor=1)
  // [NB do keep 'static' for these constexpr arrays, see issue #283]
  static constexpr fptype2 colorMatrix[ncolor][ncolor] = { { 1 } }; // 2-D array[1][1]

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
  // The fptype2 version is the default used by kernels (supporting mixed floating point mode also in blas)
  static __device__ fptype2 s_pNormalizedColorMatrix2[ncolor * ncolor];
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void createNormalizedColorMatrix()
  {
    static bool first = true;
    if( first )
    {
      first = false;
      constexpr NormalizedColorMatrix<fptype2> normalizedColorMatrix2;
      gpuMemcpyToSymbol( s_pNormalizedColorMatrix2, normalizedColorMatrix2.value, ncolor * ncolor * sizeof( fptype2 ) );
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifndef MGONGPUCPP_GPUIMPL
  void
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
  __global__ void
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

#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_HAS_NO_BLAS
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
  __global__ void
  convertD2F_Jamps( fptype2* allJampsFpt2,   // output: jamp[ncolor*2*nevt] for one specific helicity
                    const fptype* allJamps ) // input: jamp[ncolor*2*nevt] for one specific helicity
  {
    const int nevt = gridDim.x * blockDim.x;
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
    // NB! From a functional point of view, any striding will be ok here as long as ncolor*2*nevt elements are all correctly copied!
    // NB! Just in case this may be better for performance reasons, however, the same striding as in compute_jamps and cuBLAS is used here
    for( int ix2 = 0; ix2 < mgOnGpu::nx2; ix2++ )
      for( int icol = 0; icol < ncolor; icol++ )
        allJampsFpt2[ix2 * ncolor * nevt + icol * nevt + ievt] = allJamps[ix2 * ncolor * nevt + icol * nevt + ievt]; // "new1" striding
  }
#endif
#endif
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_HAS_NO_BLAS
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
  __global__ void
  convertF2D_MEs( fptype* allMEs,             // output: allMEs[nevt] for one specific helicity
                  const fptype2* allMEsFpt2 ) // input: allMEs[nevt] for one specific helicity
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
    allMEs[ievt] = allMEsFpt2[ievt];
  }
#endif
#endif
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL /* clang-format off */
#ifndef MGONGPU_HAS_NO_BLAS
  void
  color_sum_blas( fptype* allMEs,               // output: allMEs[nevt], add |M|^2 for this specific helicity
                  const fptype* allJamps,       // input: jamp[ncolor*2*nevt] for one specific helicity
                  fptype2* allBlasTmp,          // tmp: blasTmp[ncolor*2*nevt] or blasTmp[(2*ncolor*2+1)*nevt] for one specific helicity
#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
                  gpuStream_t stream,           // input: cuda stream (nullptr indicates the default stream - only used for FPTYPE=m)
#else
                  gpuStream_t /*stream*/,       // input: cuda stream (nullptr indicates the default stream - only used for FPTYPE=m)
#endif
                  gpuBlasHandle_t* pBlasHandle, // input: cuBLAS/hipBLAS handle
                  const int gpublocks,          // input: cuda gpublocks
                  const int gputhreads )        // input: cuda gputhreads
  {
    const int nevt = gpublocks * gputhreads;

    // Get the address associated with the normalized color matrix in device memory
    static fptype2* devNormColMat = nullptr;
    if( !devNormColMat ) gpuGetSymbolAddress( (void**)&devNormColMat, s_pNormalizedColorMatrix2 );

#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Mixed precision mode: need two fptype2[ncolor*2*nevt] buffers and one fptype2[nevt] buffer per helicity
    fptype2* allZtempBoth = allBlasTmp;                                  // start of first fptype2[ncolor*2*nevt] buffer
    fptype2* allJampsFpt2 = allBlasTmp + ncolor * mgOnGpu::nx2 * nevt;   // start of second fptype2[ncolor*2*nevt] buffer
    fptype2* allMEsFpt2 = allBlasTmp + 2 * ncolor * mgOnGpu::nx2 * nevt; // start of fptype2[nevt] buffer
    // Convert jamps from double to float
    gpuLaunchKernelStream( convertD2F_Jamps, gpublocks, gputhreads, stream, allJampsFpt2, allJamps );
    // Real and imaginary components
    const fptype2* allJampsReal = allJampsFpt2;
    const fptype2* allJampsImag = allJampsFpt2 + ncolor * nevt;
#else
    static_assert( std::is_same<fptype2, fptype>::value );     // sanity check
    // Standard single/double precision mode: need one fptype2[ncolor*2*nevt] buffer
    fptype2* allZtempBoth = allBlasTmp; // start of fptype2[ncolor*2*nevt] buffer
    fptype2* allMEsFpt2 = allMEs;
    // Real and imaginary components
    const fptype2* allJampsReal = allJamps;                 // this is not a cast (the two types are identical)
    const fptype2* allJampsImag = allJamps + ncolor * nevt; // this is not a cast (the two types are identical)
#endif
    // Real and imaginary components
    fptype2* allZtempReal = allZtempBoth;
    fptype2* allZtempImag = allZtempBoth + ncolor * nevt;

    // Note, new striding for cuBLAS from DeviceAccessJamp:
    // - allJamps(icol,ievt).real is allJamps[0 * ncolor * nevt + icol * nevt + ievt] // "new1"
    // - allJamps(icol,ievt).imag is allJamps[1 * ncolor * nevt + icol * nevt + ievt] // "new1"

    // Step 1: Compute Ztemp[ncolor][nevt] = ColorMatrix[ncolor][ncolor] * JampsVector[ncolor][nevt] for both real and imag
    // In this case alpha=1 and beta=0: the operation is Ztemp = alpha * ColorMatrix * JampsVector + beta * Ztemp
    fptype2 alpha1 = 1;
    fptype2 beta1 = 0;
    const int ncolorM = ncolor;
    const int nevtN = nevt;
    const int ncolorK = ncolor;
    checkGpuBlas( gpuBlasTgemm( *pBlasHandle,
                                GPUBLAS_OP_N,              // do not transpose ColMat
                                GPUBLAS_OP_T,              // transpose JampsV (new1)
                                ncolorM, nevtN, ncolorK,
                                &alpha1,
                                devNormColMat, ncolorM,    // ColMat is ncolorM x ncolorK
                                allJampsReal, nevtN,       // JampsV is nevtN x ncolorK (new1)
                                &beta1,
                                allZtempReal, ncolorM ) ); // Ztemp is ncolorM x nevtN
    checkGpuBlas( gpuBlasTgemm( *pBlasHandle,
                                GPUBLAS_OP_N,              // do not transpose ColMat
                                GPUBLAS_OP_T,              // transpose JampsV (new1)
                                ncolorM, nevtN, ncolorK,
                                &alpha1,
                                devNormColMat, ncolorM,    // ColMat is ncolorM x ncolorK
                                allJampsImag, nevtN,       // JampsV is nevtN x ncolorK (new1)
                                &beta1,
                                allZtempImag, ncolorM ) ); // Ztemp is ncolorM x nevtN

    // Step 2: For each ievt, compute the dot product of JampsVector[ncolor][ievt] dot tmp[ncolor][ievt]
    // In this case alpha=1 and beta=1: the operation is ME = alpha * ( Tmp dot JampsVector ) + beta * ME
    // Use cublasSgemmStridedBatched to perform these batched dot products in one call
    fptype2 alpha2 = 1;
    fptype2 beta2 = 1;
    checkGpuBlas( gpuBlasTgemmStridedBatched( *pBlasHandle,
                                              GPUBLAS_OP_N,                   // do not transpose JampsV (new1)
                                              GPUBLAS_OP_N,                   // do not transpose Tmp
                                              1, 1, ncolor,                   // result is 1x1 (dot product)
                                              &alpha2,
                                              allJampsReal, nevt, 1,          // allJamps is nevt x ncolor, stride 1 for each ievt column (new1)
                                              allZtempReal, ncolor, ncolor,   // allZtemp is ncolor x nevt, with stride ncolor for each ievt column
                                              &beta2,
                                              allMEsFpt2, 1, 1,               // output is a 1x1 result for each "batch" (i.e. for each ievt)
                                              nevt ) );                       // there are nevt "batches"
    checkGpuBlas( gpuBlasTgemmStridedBatched( *pBlasHandle,
                                              GPUBLAS_OP_N,                   // do not transpose JampsV (new1)
                                              GPUBLAS_OP_N,                   // do not transpose Tmp
                                              1, 1, ncolor,                   // result is 1x1 (dot product)
                                              &alpha2,
                                              allJampsImag, nevt, 1,          // allJamps is nevt x ncolor, stride 1 for each ievt column (new1)
                                              allZtempImag, ncolor, ncolor,   // allZtemp is ncolor x nevt, with stride ncolor for each ievt column
                                              &beta2,
                                              allMEsFpt2, 1, 1,               // output is a 1x1 result for each "batch" (i.e. for each ievt)
                                              nevt ) );                       // there are nevt "batches"

#if defined MGONGPU_FPTYPE_DOUBLE and defined MGONGPU_FPTYPE2_FLOAT
    // Convert MEs from float to double
    gpuLaunchKernelStream( convertF2D_MEs, gpublocks, gputhreads, stream, allMEs, allMEsFpt2 );
#endif
  }
#endif /* clang-format on */
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void
  color_sum_gpu( fptype* allMEs,               // output: allMEs[nevt], add |M|^2 for this specific helicity
                 const fptype* allJamps,       // input: jamp[ncolor*2*nevt] for one specific helicity
                 fptype2* allBlasTmp,          // tmp: blasTmp[ncolor*2*nevt] or blasTmp[(2*ncolor*2+1)*nevt] for one specific helicity
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
      color_sum_blas( allMEs, allJamps, allBlasTmp, stream, pBlasHandle, gpublocks, gputhreads );
    }
#endif
  }
#endif

  //--------------------------------------------------------------------------

} // end namespace
