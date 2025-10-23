// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Oct 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#ifndef diagrams_header_H
#define diagrams_header_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "CPPProcess.h"

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  constexpr int nw6 = CPPProcess::nw6;       // dimensions of each wavefunction (HELAS KEK 91-11): e.g. 6 for e+ e- -> mu+ mu- (fermions and vectors)
  constexpr int nwf = CPPProcess::nwf;       // #wavefunctions = #external (npar) + #internal: e.g. 5 for e+ e- -> mu+ mu- (1 internal is gamma or Z)
  constexpr int ncolor = CPPProcess::ncolor; // the number of leading colors

  using Parameters_sm_dependentCouplings::ndcoup;   // #couplings that vary event by event (depend on running alphas QC>
  using Parameters_sm_independentCouplings::nicoup; // #couplings that are fixed for all events (do not depend on runni>

#ifdef __CUDACC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // e.g. <<warning #177-D: variable "nevt" was declared but never referenced>>
#endif
  constexpr int nIPD = CPPProcess::nIPD; // SM independent parameters
  constexpr int nIPC = CPPProcess::nIPC; // SM independent couplings
#ifdef __CUDACC__
#pragma nv_diagnostic pop
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  // Encapsulate here (rather than in MemoyAccessWavefunctions.h) the wavefunction memory layout in GPU global memory
  // *** NB: Non-trivial access in GPU global memory is only used in storeWf and retrieveWf ***
  class DeviceAccessWavefunctions
  {
  public:
    static __host__ __device__ inline cxtype&
    kernelAccessIw6( fptype* buffer, const int iw6 )
    {
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      //return reinterpret_cast<cxtype*>( buffer + ievt * CPPProcess::nw6 * mgOnGpu::nx2 )[iw6]; // OLD (non coalesced?)
      const int nevt = gridDim.x * blockDim.x;
      return *( reinterpret_cast<cxtype*>( buffer + ( iw6 * nevt + ievt ) * mgOnGpu::nx2 ) ); // NEW (coalesced?)
    }
    static __host__ __device__ inline const cxtype
    kernelAccessIw6Const( const fptype* buffer, const int iw6 )
    {
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;
      //return reinterpret_cast<const cxtype*>( buffer + ievt * CPPProcess::nw6 * mgOnGpu::nx2 )[iw6]; // OLD (non coalesced?)
      const int nevt = gridDim.x * blockDim.x;
      return *( reinterpret_cast<const cxtype*>( buffer + ( iw6 * nevt + ievt ) * mgOnGpu::nx2 ) ); // NEW (coalesced?)
    }
  };
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  inline __device__ void
  retrieveWf( const fptype* allWfs,
              cxtype w_cx[][CPPProcess::nw6],
              int nevt,
              int iwf )
  {
    using WG_ACCESS = DeviceAccessWavefunctions; // non-trivial access in global memory
    const fptype* allWfs_iwf = allWfs + iwf * nevt * CPPProcess::nw6 * mgOnGpu::nx2;
    // NB copy all 6 components (only the last 4 are used to compute amplitudes, but all 6 are needed to compute other wavefunctions)
    for( int iw6 = 0; iw6 < CPPProcess::nw6; iw6++ )
      w_cx[iwf][iw6] = WG_ACCESS::kernelAccessIw6Const( allWfs_iwf, iw6 );
  }
#endif

  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  inline __device__ void
  storeWf( fptype* allWfs,
           const cxtype w_cx[][CPPProcess::nw6],
           int nevt,
           int iwf )
  {
    using WG_ACCESS = DeviceAccessWavefunctions; // non-trivial access in global memory
    fptype* allWfs_iwf = allWfs + iwf * nevt * CPPProcess::nw6 * mgOnGpu::nx2;
    // NB copy all 6 components (only the last 4 are used to compute amplitudes, but all 6 are needed to compute other wavefunctions)
    for( int iw6 = 0; iw6 < CPPProcess::nw6; iw6++ )
      WG_ACCESS::kernelAccessIw6( allWfs_iwf, iw6 ) = w_cx[iwf][iw6];
  }
#endif

  //--------------------------------------------------------------------------
}

#endif // diagrams_header_H
