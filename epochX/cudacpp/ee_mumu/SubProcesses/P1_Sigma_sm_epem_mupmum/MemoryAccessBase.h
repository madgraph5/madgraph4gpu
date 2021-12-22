#ifndef MemoryAccessBase_H
#define MemoryAccessBase_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

// A collection of boilerplate templated functions for MemoryAccess classes
namespace MemoryAccessBase
{

  // (Non-const memory access to field from ievent)
  template<class T>
  __host__ __device__ inline
  fptype& ieventAccessField( fptype* buffer,
                             const int ievt,
                             const int ip4,
                             const int ipar )
  {
    // NB all KernelLauncher classes assume that memory access can be decomposed in this way
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( T::ieventAccessRecord( buffer, ievt ), ip4, ipar );
  }

}

#endif // MemoryAccessBase_H
