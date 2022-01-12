#ifndef MemoryAccessHelpers_H
#define MemoryAccessHelpers_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuFptypes.h"

//----------------------------------------------------------------------------

// A templated helper class that includes the boilerplate code for MemoryAccess classes
template<class T>
class MemoryAccessHelper
{
public:

  //--------------------------------------------------------------------------

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // Non-const memory access to field from ievent
  // [Signature ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecord = T::ieventAccessRecord;

  //--------------------------------------------------------------------------

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // Const memory access to event record from ievent
  // [Signature ===> fptype* ieventAccessRecordConst( const fptype* buffer, const int ievt ) <===]
  static
  __host__ __device__ inline
  fptype* ieventAccessRecordConst( const fptype* buffer,
                                   const int ievt )
  {
    return ieventAccessRecord( const_cast<fptype*>( buffer ), ievt );
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // Non-const memory access to field in an event record
  // [Signature ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
  static constexpr auto decodeRecord = T::decodeRecord;

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // Const memory access to field in an event record
  // [Signature ===> const fptype& decodeRecordConst( fptype* buffer, Ts... args ) <===]
  template<class... Ts>
  static
  __host__ __device__ inline
  const fptype& decodeRecordConst( fptype* buffer,
                                   Ts... args ) // variadic template
  {
    return T::decodeRecord( const_cast<fptype*>( buffer ), args... );
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // Non-const memory access to field from ievent
  // [Signature ===> fptype& ieventAccessField( fptype* buffer, const ievt, Ts... args ) <===]
  template<class... Ts>
  static
  __host__ __device__ inline
  fptype& ieventAccessField( fptype* buffer,
                             const int ievt,
                             Ts... args ) // variadic template
  {
    // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( T::ieventAccessRecord( buffer, ievt ), args... );
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // Const memory access to field from ievent
  // [Signature ===> const fptype& ieventAccessFieldConst( const fptype* buffer, const ievt, Ts... args ) <===]
  template<class... Ts>
  static
  __host__ __device__ inline
  const fptype& ieventAccessFieldConst( const fptype* buffer,
                                        const int ievt,
                                        Ts... args ) // variadic template
  {
    return ieventAccessField( const_cast<fptype*>( buffer ), ievt, args... );
  }

};

//----------------------------------------------------------------------------

// A templated helper class that includes the boilerplate code for KernelAccess classes
template<class T, bool onDevice>
class KernelAccessHelper : public MemoryAccessHelper<T>
{
public:

  //--------------------------------------------------------------------------

  // Locate an event record (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal)
  // Non-const memory access to event record from kernel
  // [Signature ===> fptype* kernelAccessRecord( fptype* buffer ) <===]
  static
  __host__ __device__ inline
  fptype* kernelAccessRecord( fptype* buffer )
  {
    //if constexpr ( !onDevice ) // FIXME! enable this when we move to nvcc supporting c++17
    if ( !onDevice )
    {
      return T::ieventAccessRecord( buffer, 0 );
    }
    else
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "kernelAccessRecord: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      return T::ieventAccessRecord( buffer, ievt ); // NB fptype and fptype_sv coincide for CUDA
#else
      throw std::runtime_error( "kernelAccessRecord on device is only implemented in CUDA" );
#endif
    }
  }

  //--------------------------------------------------------------------------
 
  // Locate an event record (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal)
  // Const memory access to event record from kernel
  // [Signature ===> fptype* kernelAccessRecordConst( const fptype* buffer ) <===]
  static
  __host__ __device__ inline
  fptype* kernelAccessRecordConst( const fptype* buffer )
  {
    return kernelAccessRecord( const_cast<fptype*>( buffer ) );
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // Non-const memory access to field from kernel
  // [Signature ===> fptype& kernelAccessField( fptype* buffer, Ts... args ) <===]
  template<class... Ts>
  static
  __host__ __device__ inline
  fptype& kernelAccessField( fptype* buffer,
                             Ts... args ) // variadic template
  {
    // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( kernelAccessRecord( buffer ), args... );
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // Const memory access to field from kernel
  // [Signature ===> const fptype& kernelAccessFieldConst( const fptype* buffer, Ts... args ) <===]
  template<class... Ts>
  static
  __host__ __device__ inline
  const fptype& kernelAccessFieldConst( const fptype* buffer,
                                        Ts... args ) // variadic template
  {
    return kernelAccessField( const_cast<fptype*>( buffer ), args... );
  }

  //--------------------------------------------------------------------------

};

#endif // MemoryAccessHelpers_H
