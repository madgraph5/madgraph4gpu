#ifndef MemoryAccessHelpers_H
#define MemoryAccessHelpers_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

//----------------------------------------------------------------------------

// A templated helper class that includes the boilerplate code for MemoryAccess classes
template<class T>
class MemoryAccessHelper
{
public:

  //--------------------------------------------------------------------------

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccessRecord = T::ieventAccessRecord;

  //--------------------------------------------------------------------------

  // (Const memory access to event record from ievent)
  static
  __host__ __device__ inline
  fptype* ieventAccessRecordConst( const fptype* buffer,
                                   const int ievt )
  {
    return ieventAccessRecord( const_cast<fptype*>( buffer ), ievt );
  }

  //--------------------------------------------------------------------------

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecord = T::decodeRecord;

  //--------------------------------------------------------------------------

  // (Const memory access to field in an event record)
  template<class... Ts> // variadic template
  static
  __host__ __device__ inline
  const fptype& decodeRecordConst( fptype* buffer,
                                   Ts... args )
  {
    return T::decodeRecord( const_cast<fptype*>( buffer ), args... );
  }

  //--------------------------------------------------------------------------

  // (Non-const memory access to field from ievent)
  template<class... Ts> // variadic template
  static
  __host__ __device__ inline
  fptype& ieventAccessField( fptype* buffer,
                             const int ievt,
                             Ts... args )
  {
    // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( T::ieventAccessRecord( buffer, ievt ), args... );
  }

  //--------------------------------------------------------------------------

  // (Const memory access to field from ievent)
  template<class... Ts> // variadic template
  static
  __host__ __device__ inline
  const fptype& ieventAccessFieldConst( const fptype* buffer,
                                        const int ievt,
                                        Ts... args )
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

  // Locate an event record (output) in a memory buffer (input) from an implicit event-indexing mechanism in the kernel
  // (Non-const memory access to event record from kernel)
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

  // (Const memory access to event record from kernel)
  static
  __host__ __device__ inline
  fptype* kernelAccessRecordConst( const fptype* buffer )
  {
    return kernelAccessRecord( const_cast<fptype*>( buffer ) );
  }

  //--------------------------------------------------------------------------

  // (Non-const memory access to field from kernel)
  template<class... Ts> // variadic template
  static
  __host__ __device__ inline
  fptype& kernelAccessField( fptype* buffer,
                             Ts... args )
  {
    // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( kernelAccessRecord( buffer ), args... );
  }

  //--------------------------------------------------------------------------

  // (Const memory access to field from kernel)
  template<class... Ts> // variadic template
  static
  __host__ __device__ inline
  const fptype& kernelAccessFieldConst( const fptype* buffer,
                                        Ts... args )
  {
    return kernelAccessField( const_cast<fptype*>( buffer ), args... );
  }

  //--------------------------------------------------------------------------

};

#endif // MemoryAccessHelpers_H
