// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Oct 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#ifndef diagrams_header_H
#define diagrams_header_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "CPPProcess.h"
#include "MemoryAccessWavefunctions.h"

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
{
  //--------------------------------------------------------------------------

  inline __device__ void
  retrieveWf( const fptype* allWfs,
              cxtype w_cx[][CPPProcess::nw6],
              int nevt,
              int iwf )
  {
    using WG_ACCESS = DeviceAccessWavefunctions; // non-trivial access in global memory
    const fptype* allWfs_iwf = allWfs + iwf * nevt * CPPProcess::nw6 * mgOnGpu::nx2;
    const cxtype* wcx_iwf = WG_ACCESS::kernelAccessConst( allWfs_iwf );
    for( int iw6 = 0; iw6 < CPPProcess::nw6; iw6++ ) // FIXME: only need 4 out of 6?
      w_cx[iwf][iw6] = wcx_iwf[iw6];
  }

  //--------------------------------------------------------------------------

  inline __device__ void
  storeWf( fptype* allWfs,
           const cxtype w_cx[][CPPProcess::nw6],
           int nevt,
           int iwf )
  {
    using WG_ACCESS = DeviceAccessWavefunctions; // non-trivial access in global memory
    fptype* allWfs_iwf = allWfs + iwf * nevt * CPPProcess::nw6 * mgOnGpu::nx2;
    cxtype* wcx_iwf = WG_ACCESS::kernelAccess( allWfs_iwf );
    for( int iw6 = 0; iw6 < CPPProcess::nw6; iw6++ ) // FIXME: only need 4 out of 6?
      wcx_iwf[iw6] = w_cx[iwf][iw6];
  }

  //--------------------------------------------------------------------------
}
#endif

#endif // diagrams_header_H
