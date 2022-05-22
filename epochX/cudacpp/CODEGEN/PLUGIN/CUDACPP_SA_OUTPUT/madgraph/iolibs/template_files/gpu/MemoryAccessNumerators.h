#ifndef MemoryAccessNumerators_H
#define MemoryAccessNumerators_H 1

#include "MemoryAccessGs.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for numerators
// This implementation reuses the plain ARRAY[nevt] implementation of MemoryAccessGs

typedef KernelAccessGs<false> HostAccessNumerators;
typedef KernelAccessGs<true> DeviceAccessNumerators;

//----------------------------------------------------------------------------

#endif // MemoryAccessNumerators_H
