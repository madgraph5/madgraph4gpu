// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (May 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef MemoryAccessDenominators_H
#define MemoryAccessDenominators_H 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL

#include "MemoryAccessGs.h"

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //----------------------------------------------------------------------------

  // A class describing the internal layout of memory buffers for denominators
  // This implementation reuses the plain ARRAY[nevt] implementation of MemoryAccessGs

  typedef KernelAccessGs<false> HostAccessDenominators;
  typedef KernelAccessGs<true> DeviceAccessDenominators;

  //----------------------------------------------------------------------------

} // end namespace mg5amcGpu/mg5amcCpu

#endif
#endif // MemoryAccessDenominators_H
