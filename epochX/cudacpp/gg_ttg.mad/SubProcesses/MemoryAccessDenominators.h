#ifndef MemoryAccessDenominators_H
#define MemoryAccessDenominators_H 1

#include "MemoryAccessGs.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for denominators
// This implementation reuses the plain ARRAY[nevt] implementation of MemoryAccessGs

typedef KernelAccessGs<false> HostAccessDenominators;
typedef KernelAccessGs<true> DeviceAccessDenominators;

//----------------------------------------------------------------------------

#endif // MemoryAccessDenominators_H
