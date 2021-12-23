#ifndef MemoryAccessBase_H
#define MemoryAccessBase_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

//----------------------------------------------------------------------------

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
  template<class... Ts>
  static
  __host__ __device__ inline
  fptype& ieventAccessFIELD( fptype* buffer,
                             const int ievt,
                             Ts... args )
  {
    // NB all KernelLauncher classes assume that memory access can be decomposed in this way
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( T::ieventAccessRecord( buffer, ievt ), args... );
  }

  /*
  static constexpr auto ieventAccessFIELD2 = ieventAccessFIELD<const int, const int>; // builds or not builds depending on below

  static
  __host__ __device__ inline
  fptype& ieventAccessField( fptype* buffer,
                             const int ievt,
                             const int ip4,
                             const int ipar )
  {
    //return ieventAccessFIELD<const int, const int>( buffer, ievt, ip4, ipar ); // builds
    return ieventAccessFIELD2( buffer, ievt, ip4, ipar ); // triggers above line not to build
    //return *buffer; // builds
  }
  */

  //--------------------------------------------------------------------------

  /*
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
  */

};

//----------------------------------------------------------------------------

// A templated base class that includes the boilerplate code for KernelAccess classes
template<class T, bool onDevice>
class KernelAccessBase : public MemoryAccessBase<T>
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
  fptype* kernelAccessConstRecord( const fptype* buffer )
  {
    return kernelAccessRecord( const_cast<fptype*>( buffer ) );
  }

  //--------------------------------------------------------------------------

  // (Non-const memory access to field from kernel)
  static
  __host__ __device__ inline
  fptype& kernelAccessField( fptype* buffer,
                             const int ip4,
                             const int ipar )
  {
    // NB all KernelLauncher classes assume that memory access can be decomposed in this way
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return T::decodeRecord( kernelAccessRecord( buffer ), ip4, ipar );
  }

  //--------------------------------------------------------------------------

  // (Const memory access to field from ievent)
  static
  __host__ __device__ inline
  const fptype& kernelAccessConstField( const fptype* buffer,
                                        const int ip4,
                                        const int ipar )
  {
    return kernelAccessField( const_cast<fptype*>( buffer ), ip4, ipar );
  }

  //--------------------------------------------------------------------------

};

#endif // MemoryAccessBase_H
