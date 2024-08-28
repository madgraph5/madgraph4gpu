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

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  __device__ void
  linker_VVV1P0_1( const fptype allV2[],
                   const fptype allV3[],
                   const fptype allCOUP[],
                   const double Ccoeff,
                   const fptype M1,
                   const fptype W1,
                   fptype allV1[] )
  {
    return VVV1P0_1<W_ACCESS, CD_ACCESS>( allV2, allV3, allCOUP, Ccoeff, M1, W1, allV1 );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], S2[6], S3[6]
  __device__ void
  linker_VSS1_0( const fptype allV1[],
                 const fptype allS2[],
                 const fptype allS3[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 fptype allvertexes[] )
  {
    return VSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( allV1, allS2, allS3, allCOUP, Ccoeff, allvertexes );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'S2[6]' from the input wavefunctions V1[6], S3[6]
  __device__ void
  linker_VSS1_2( const fptype allV1[],
                 const fptype allS3[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 const fptype M2,
                 const fptype W2,
                 fptype allS2[] )
  {
    return VSS1_2<W_ACCESS, CD_ACCESS>( allV1, allS3, allCOUP, Ccoeff, M2, W2, allS2 );
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'S3[6]' from the input wavefunctions V1[6], S2[6]
  __device__ void
  linker_VSS1_3( const fptype allV1[],
                 const fptype allS2[],
                 const fptype allCOUP[],
                 const double Ccoeff,
                 const fptype M3,
                 const fptype W3,
                 fptype allS3[] )
  {
    return VSS1_3<W_ACCESS, CD_ACCESS>( allV1, allS2, allCOUP, Ccoeff, M3, W3, allS3 );
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions V1[6], V2[6], S3[6], S4[6]
  __device__ void
  linker_VVSS1_0( const fptype allV1[],
                  const fptype allV2[],
                  const fptype allS3[],
                  const fptype allS4[],
                  const fptype allCOUP[],
                  const double Ccoeff,
                  fptype allvertexes[] )
  {
    return VVSS1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( allV1, allV2, allS3, allS4, allCOUP, Ccoeff, allvertexes );
  }

  //--------------------------------------------------------------------------
}
#endif
