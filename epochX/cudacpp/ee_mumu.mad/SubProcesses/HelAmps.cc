// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Aug 2024) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2024) for the MG5aMC CUDACPP plugin.

#ifdef MGONGPU_LINKER_HELAMPS

#include "HelAmps_sm.h"

// -----------------------------------------------------------------------------
// *** NB: this implementation class depends on MemoryAccessMomenta,
// *** where the AOSOA definition depends on CPPProcess::npar,
// *** which may be different in different P* subprocess directories:
// *** therefore this class is presently hosted and compiled in each P*
// -----------------------------------------------------------------------------

#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessGs.h"
#include "MemoryAccessMatrixElements.h"
#include "MemoryAccessMomenta.h"
#include "MemoryAccessWavefunctions.h"

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "MemoryAccessDenominators.h"
#include "MemoryAccessNumerators.h"
#endif

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

#ifdef MGONGPUCPP_GPUIMPL
  using M_ACCESS = DeviceAccessMomenta;         // non-trivial access: buffer includes all events
  using E_ACCESS = DeviceAccessMatrixElements;  // non-trivial access: buffer includes all events
  using W_ACCESS = DeviceAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
  using A_ACCESS = DeviceAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
  using CD_ACCESS = DeviceAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
  using CI_ACCESS = DeviceAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  using NUM_ACCESS = DeviceAccessNumerators;   // non-trivial access: buffer includes all events
  using DEN_ACCESS = DeviceAccessDenominators; // non-trivial access: buffer includes all events
#endif
#else
  using namespace ::mg5amcCpu;
  using M_ACCESS = HostAccessMomenta;         // non-trivial access: buffer includes all events
  using E_ACCESS = HostAccessMatrixElements;  // non-trivial access: buffer includes all events
  using W_ACCESS = HostAccessWavefunctions;   // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
  using A_ACCESS = HostAccessAmplitudes;      // TRIVIAL ACCESS (no kernel splitting yet): buffer for one event
  using CD_ACCESS = HostAccessCouplings;      // non-trivial access (dependent couplings): buffer includes all events
  using CI_ACCESS = HostAccessCouplingsFixed; // TRIVIAL access (independent couplings): buffer for one event
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
  using NUM_ACCESS = HostAccessNumerators;    // non-trivial access: buffer includes all events
  using DEN_ACCESS = HostAccessDenominators;  // non-trivial access: buffer includes all events
#endif
#endif

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ void
  linker_FFV1_0( const fptype allF1[],
                 const fptype allF2[],
                 const fptype allV3[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 fptype allvertexes[] )
  {
    return FFV1_0<W_ACCESS, A_ACCESS, CI_ACCESS>( allF1, allF2, allV3, allCOUP, Ccoeff, allvertexes );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ void
  linker_FFV1P0_3( const fptype allF1[],
                   const fptype allF2[],
                   const fptype allCOUP[],
                   const double Ccoeff,
                   const fptype M3,
                   const fptype W3,
                   fptype allV3[] )
  {
    return FFV1P0_3<W_ACCESS, CI_ACCESS>( allF1, allF2, allCOUP, Ccoeff, M3, W3, allV3 );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ void
  linker_FFV2_0( const fptype allF1[],
                 const fptype allF2[],
                 const fptype allV3[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 fptype allvertexes[] )
  {
    return FFV2_0<W_ACCESS, A_ACCESS, CI_ACCESS>( allF1, allF2, allV3, allCOUP, Ccoeff, allvertexes );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ void
  linker_FFV2_3( const fptype allF1[],
                 const fptype allF2[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 const fptype M3,
                 const fptype W3,
                 fptype allV3[] )
  {
    return FFV2_3<W_ACCESS, CI_ACCESS>( allF1, allF2, allCOUP, Ccoeff, M3, W3, allV3 );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ void
  linker_FFV4_0( const fptype allF1[],
                 const fptype allF2[],
                 const fptype allV3[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 fptype allvertexes[] )
  {
    return FFV4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( allF1, allF2, allV3, allCOUP, Ccoeff, allvertexes );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ void
  linker_FFV4_3( const fptype allF1[],
                 const fptype allF2[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 const fptype M3,
                 const fptype W3,
                 fptype allV3[] )
  {
    return FFV4_3<W_ACCESS, CI_ACCESS>( allF1, allF2, allCOUP, Ccoeff, M3, W3, allV3 );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__ void
  linker_FFV2_4_0( const fptype allF1[],
                   const fptype allF2[],
                   const fptype allV3[],
                   const fptype allCOUP[],
                   const double Ccoeff,
                   fptype allvertexes[] )
  {
    return FFV2_4_0<W_ACCESS, A_ACCESS, CI_ACCESS>( allF1, allF2, allV3, allCOUP, Ccoeff, allvertexes );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__ void
  linker_FFV2_4_3( const fptype allF1[],
                   const fptype allF2[],
                   const fptype allCOUP[],
                   const double Ccoeff,
                   const fptype M3,
                   const fptype W3,
                   fptype allV3[] )
  {
    return FFV2_4_3<W_ACCESS, CI_ACCESS>( allF1, allF2, allCOUP, Ccoeff, M3, W3, allV3 );
  }

  //--------------------------------------------------------------------------
}
#endif
