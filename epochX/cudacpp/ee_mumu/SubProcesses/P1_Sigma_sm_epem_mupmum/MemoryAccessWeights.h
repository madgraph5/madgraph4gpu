#ifndef MemoryAccessWeights_H
#define MemoryAccessWeights_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for weights
// This implementation uses a plain ARRAY[nevt]
// [If many implementations are used, a suffix _ARRAYv1 should be appended to the class name]
class MemoryAccessWeights//_ARRAYv1
{
public:

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Indexed access (WITH an explicit event number) to the memory buffer for weights
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific weight for one event, given its event number
  // (Non-const memory access)
  static
  __host__ __device__ inline
  fptype& ieventAccess( fptype* buffer,
                        const int ievt )
  {
    return buffer[ievt]; // ARRAY[nevt]
  }

  // (Const memory access)
  static
  __host__ __device__ inline
  const fptype& ieventConstAccess( fptype* buffer,
                                   const int ievt )
  {
    return ieventAccess( const_cast<fptype*>( buffer ), ievt );
  }

};

//----------------------------------------------------------------------------

// A class describing kernel access to memory buffers for weights on a CPU host or on a GPU device
template<bool onDevice>
class KernelAccessWeights: public MemoryAccessWeights
{
public:

  // =========================================================================
  // *** Pattern: kernelAccessInd1..IndN( buffer [, ind1... indN] )        ***
  // =========================================================================

  // Kernel access (WITHOUT an explicit event number) to the memory buffer for weights
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific weight for one event, given its event number
  // (Non-const memory access)
  static
  __host__ __device__ inline
  fptype& kernelAccess( fptype* buffer )
  {
    //if constexpr ( !onDevice ) // FIXME! enable this when we move to nvcc supporting c++17
    if ( !onDevice )
    {
      return ieventAccess( buffer, 0 );
    }
    else
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "kernelAccessIp4Ipar: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      return ieventAccess( buffer, ievt ); // NB fptype and fptype_sv coincide for CUDA
#else
      throw std::runtime_error( "KernelAccessWeights on device is only implemented in CUDA" );
#endif
    }
  }

  // (Const memory access)
  static
  __host__ __device__ inline
  const fptype& kernelConstAccess( fptype* buffer )
  {
    return kernelAccess( const_cast<fptype*>( buffer ) );
  }

};

//----------------------------------------------------------------------------

typedef KernelAccessWeights<false> HostAccessWeights;
typedef KernelAccessWeights<true> DeviceAccessWeights;

//----------------------------------------------------------------------------

#endif // MemoryAccessWeights_H
