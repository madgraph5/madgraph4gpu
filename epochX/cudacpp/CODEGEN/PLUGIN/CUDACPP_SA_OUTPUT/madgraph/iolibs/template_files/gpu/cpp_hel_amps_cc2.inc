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
%(function_definitions2)s}
#endif
