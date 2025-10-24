// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#ifndef COLOR_SUM_H
#define COLOR_SUM_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuVectors.h"

#include "CPPProcess.h"
#include "GpuAbstraction.h"

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  class DeviceAccessJamp
  {
  public:
    static __device__ inline cxtype_ref
    kernelAccessIcol( fptype* buffer, const int icol )
    {
      const int ncolor = CPPProcess::ncolor; // the number of leading colors
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      // Original "old" striding for CUDA kernels: ncolor separate 2*nevt matrices for each color (ievt last)
      //return cxtype_ref( buffer[icol * 2 * nevt + ievt], buffer[icol * 2 * nevt + nevt + ievt] ); // "old"
      // New "new1" striding for cuBLAS: two separate ncolor*nevt matrices for each of real and imag (ievt last)
      // The "new1" striding is now used for both HASBLAS=hasBlas and hasNoBlas builds and for both CUDA kernels and cuBLAS
      return cxtype_ref( buffer[0 * ncolor * nevt + icol * nevt + ievt], buffer[1 * ncolor * nevt + icol * nevt + ievt] ); // "new1"
      // New "new2" striding for cuBLAS: two separate ncolor*nevt matrices for each of real and imag (icol last)
      //return cxtype_ref( buffer[0 * nevt * ncolor + ievt * ncolor + icol], buffer[1 * nevt * ncolor + ievt * ncolor + icol] ); // "new2"
    }
    static __device__ inline const cxtype
    kernelAccessIcolConst( const fptype* buffer, const int icol )
    {
      const int ncolor = CPPProcess::ncolor; // the number of leading colors
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      // Original "old" striding for CUDA kernels: ncolor separate 2*nevt matrices for each color (ievt last)
      //return cxtype( buffer[icol * 2 * nevt + ievt], buffer[icol * 2 * nevt + nevt + ievt] ); // "old"
      // New "new1" striding for cuBLAS: two separate ncolor*nevt matrices for each of real and imag (ievt last)
      // The "new1" striding is now used for both HASBLAS=hasBlas and hasNoBlas builds and for both CUDA kernels and cuBLAS
      return cxtype( buffer[0 * ncolor * nevt + icol * nevt + ievt], buffer[1 * ncolor * nevt + icol * nevt + ievt] ); // "new1"
      // New "new2" striding for cuBLAS: two separate ncolor*nevt matrices for each of real and imag (icol last)
      //return cxtype( buffer[0 * nevt * ncolor + ievt * ncolor + icol], buffer[1 * nevt * ncolor + ievt * ncolor + icol] ); // "new2"
    }
  };
#else
  class HostAccessJamp
  {
  public:
    static inline cxtype_sv&
    kernelAccessIcol( cxtype_sv* buffer, const int icol )
    {
      return buffer[icol];
    }
    static inline cxtype_sv&
    kernelAccessIcol( fptype* buffer, const int icol )
    {
      return reinterpret_cast<cxtype_sv*>( buffer )[icol];
    }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  void createNormalizedColorMatrix();
#endif

  //--------------------------------------------------------------------------

#ifndef MGONGPUCPP_GPUIMPL
  void
  color_sum_cpu( fptype* allMEs,              // output: allMEs[nevt], add |M|^2 for this specific helicity
                 const cxtype_sv* allJamp_sv, // input: jamp_sv[ncolor] (float/double) or jamp_sv[2*ncolor] (mixed) for one specific helicity
                 const int ievt0 );           // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
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
                 const int gputhreads );       // input: cuda gputhreads
#endif

  //--------------------------------------------------------------------------
}

#endif // COLOR_SUM_H
