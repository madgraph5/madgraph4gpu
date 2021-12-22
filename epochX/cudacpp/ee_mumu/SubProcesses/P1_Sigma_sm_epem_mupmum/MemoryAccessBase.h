#ifndef MemoryAccessBase_H
#define MemoryAccessBase_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

// A templated base class that includes the boilerplate code for MemoryAccess classes
template<class T>
class MemoryAccessBase
{
public:

  //--------------------------------------------------------------------------

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccessRecord = T::ieventAccessRecord;
  
  //--------------------------------------------------------------------------

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecord = T::decodeRecord;

  //--------------------------------------------------------------------------

  // (Const memory access to event record from ievent)
  static
  __host__ __device__ inline
  fptype* ieventAccessConstRecord( const fptype* buffer,
                                   const int ievt )
  {
    return ieventAccessRecord( const_cast<fptype*>( buffer ), ievt );
  }

  //--------------------------------------------------------------------------

  // (Non-const memory access to field from ievent)
  static
  __host__ __device__ inline
  fptype& ieventAccessField( fptype* buffer,
                             const int ievt,
                             const int ip4,
                             const int ipar )
  {
    // NB all KernelLauncher classes assume that memory access can be decomposed in this way
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return decodeRecord( ieventAccessRecord( buffer, ievt ), ip4, ipar );
  }

  //--------------------------------------------------------------------------

  // (Const memory access to field from ievent)
  static
  __host__ __device__ inline
  const fptype& ieventAccessConstField( const fptype* buffer,
                                        const int ievt,
                                        const int ip4,
                                        const int ipar )
  {
    return ieventAccessField( const_cast<fptype*>( buffer ), ievt, ip4, ipar );
  }

  //--------------------------------------------------------------------------

};

#endif // MemoryAccessBase_H
