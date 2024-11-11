// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef MemoryAccessWavefunctions_H
#define MemoryAccessWavefunctions_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "CPPProcess.h"

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //----------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  class DeviceAccessWavefunctions
  {
  public:
    static __host__ __device__ inline cxtype_sv*
    kernelAccess( fptype* buffer )
    {
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return reinterpret_cast<cxtype*>( buffer + ievt * CPPProcess::nw6 * mgOnGpu::nx2 );
    }
    static __host__ __device__ inline const cxtype_sv*
    kernelAccessConst( const fptype* buffer )
    {
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      return reinterpret_cast<const cxtype*>( buffer + ievt * CPPProcess::nw6 * mgOnGpu::nx2 );
    }
  };
#endif
  
  //----------------------------------------------------------------------------

  class HostAccessWavefunctions
  {
  public:
    static __host__ __device__ inline cxtype_sv*
    kernelAccess( fptype* buffer )
    {
      return reinterpret_cast<cxtype_sv*>( buffer );
    }
    static __host__ __device__ inline const cxtype_sv*
    kernelAccessConst( const fptype* buffer )
    {
      return reinterpret_cast<const cxtype_sv*>( buffer );
    }
  };

  //----------------------------------------------------------------------------

} // end namespace mg5amcGpu/mg5amcCpu

#endif // MemoryAccessWavefunctions_H
