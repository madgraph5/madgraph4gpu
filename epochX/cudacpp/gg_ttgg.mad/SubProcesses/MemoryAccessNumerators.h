// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (May 2022) for the MG5aMC CUDACPP plugin.

#ifndef MemoryAccessNumerators_H
#define MemoryAccessNumerators_H 1
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL

#include "MemoryAccessGs.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for numerators
// This implementation reuses the plain ARRAY[nevt] implementation of MemoryAccessGs

typedef KernelAccessGs<false> HostAccessNumerators;
typedef KernelAccessGs<true> DeviceAccessNumerators;

//----------------------------------------------------------------------------

#endif
#endif // MemoryAccessNumerators_H
