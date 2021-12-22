#ifndef MemoryAccessRandomNumbers_H
#define MemoryAccessRandomNumbers_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

// A class describing the internal layout of memory buffers for random numbers
// This implementation uses an AOSOA[npagR][nparf][np4][neppR] where nevt=npagR*neppR
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessRandomNumbers//_AOSOAv1
{
public:

  // Number of Events Per Page in the random number AOSOA memory layout
  // *** NB Different values of neppR lead to different physics results: the ***
  // *** same 1d array is generated, but it is interpreted in different ways ***
  static constexpr int neppR = 8; // HARDCODED TO GIVE ALWAYS THE SAME PHYSICS RESULTS!
  //static constexpr int neppR = 1; // AOS (tests of sectors/requests)

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Indexed access (WITH an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the random number for a specific 4-momenta component for a specific particle in one event, given its event number
  // (Const memory access)
  static
  __device__ inline
  const fptype& ieventConstAccessIp4Iparf( const fptype* buffer,
                                           const int ievt,
                                           const int ip4,
                                           const int iparf )
  {
    using mgOnGpu::np4;
    using mgOnGpu::nparf;
    const int ipagR = ievt/neppR; // #event "R-page"
    const int ieppR = ievt%neppR; // #event in the current event R-page
    //printf( "%2d %2d %8d %8.3f\n", iparf, 0, ievt, buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR] );
    return buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR]; // AOSOA[ipagR][iparf][ip4][ieppR]
  }

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Kernel access (WITHOUT an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the random number for a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  static
  __device__ inline
  const fptype& kernelConstAccessIp4Iparf( const fptype* buffer,
                                           const int ip4,
                                           const int iparf )
  {
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "kernelCoonstAccessIp4Iparf: ievt=%d threadId=%d\n", ievt, threadIdx.x );
    return ieventConstAccessIp4Iparf( buffer, ievt, ip4, iparf ); // NB fptype and fptype_sv coincide for CUDA
#else
    return ieventConstAccessIp4Iparf( buffer, 0, ip4, iparf );
#endif
  }

  //--------------------------------------------------------------------------

};

#endif // MemoryAccessRandomNumbers_H
