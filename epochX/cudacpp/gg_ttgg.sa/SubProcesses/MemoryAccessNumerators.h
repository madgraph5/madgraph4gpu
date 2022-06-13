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
