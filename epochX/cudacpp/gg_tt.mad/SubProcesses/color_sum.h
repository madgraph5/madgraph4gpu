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
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return cxtype_ref( buffer[icol * 2 * nevt + ievt], buffer[icol * 2 * nevt + nevt + ievt] );
    }
    static __device__ inline const cxtype
    kernelAccessIcolConst( const fptype* buffer, const int icol )
    {
      const int nevt = gridDim.x * blockDim.x;
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return cxtype( buffer[icol * 2 * nevt + ievt], buffer[icol * 2 * nevt + nevt + ievt] );
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

  __global__ void /* clang-format off */
  color_sum( fptype* allMEs,              // output: allMEs[nevt], add |M|^2 for this specific helicity
#ifdef MGONGPUCPP_GPUIMPL
             const fptype* allJamps       // input: jamp[ncolor*2*nevt] for one specific helicity
#else
             const cxtype_sv* jamp_sv,    // input: jamp_sv[ncolor] (f/d) or [2*ncolor] (m) for SIMD event page(s) ievt00 and helicity ihel
             const int ievt00             // input: first event number in current C++ event page (for CUDA, ievt depends on threadid)
#endif
             ); /* clang-format on */

  //--------------------------------------------------------------------------
}

#endif // COLOR_SUM_H
