#ifndef MemoryAccessRandomNumbers_H
#define MemoryAccessRandomNumbers_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

// A class describing the internal layout of memory buffers for random numbers
// This implementation uses an AOSOA[npagR][nparf][np4][neppR] where nevt=npagR*neppR
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
//class MemoryAccessRandomNumbers//_AOSOAv1
namespace MemoryAccessRandomNumbers//_AOSOAv1
{
  //public:

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Indexed access (WITH an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  //static
  __device__ inline
  fptype& ieventAccessIp4Iparf( fptype* buffer,
                                const int ievt,
                                const int ip4,
                                const int iparf )
  {
    using mgOnGpu::np4;
    using mgOnGpu::nparf;
    constexpr int neppR = mgOnGpu::neppR; // AOSOA layout: constant at compile-time
    const int ipagR = ievt/neppR; // #event "R-page"
    const int ieppR = ievt%neppR; // #event in the current event R-page
    //printf( "%2d %2d %8d %8.3f\n", iparf, 0, ievt, buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR] );
    return buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR]; // AOSOA[ipagR][iparf][ip4][ieppR]
  }

  // (Const memory access)
  //static
  __device__ inline
  const fptype& ieventConstAccessIp4Iparf( const fptype* buffer,
                                           const int ievt,
                                           const int ip4,
                                           const int iparf )
  {
    return ieventAccessIp4Iparf( const_cast<fptype*>( buffer ), ievt, ip4, iparf );
  }

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Kernel access (WITHOUT an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  //static
  __device__ inline
  fptype& kernelAccessIp4Iparf( fptype* buffer,
                                const int ip4,
                                const int iparf )
  {
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "kernelAccessIp4Iparf: ievt=%d threadId=%d\n", ievt, threadIdx.x );
    return ieventAccessIp4Iparf( buffer, ievt, ip4, iparf ); // NB fptype and fptype_sv coincide for CUDA
#else
    return ieventAccessIp4Iparf( buffer, 0, ip4, iparf );
#endif
  }

  // (Const memory access)
  //static
  __device__ inline
  const fptype& kernelConstAccessIp4Iparf( const fptype* buffer,
                                           const int ip4,
                                           const int iparf )
  {
    return kernelAccessIp4Iparf( const_cast<fptype*>( buffer ), ip4, iparf );
  }

  //--------------------------------------------------------------------------

};

#endif // MemoryAccessRandomNumbers_H
