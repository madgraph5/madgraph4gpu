#ifndef MemoryAccessMomenta_H
#define MemoryAccessMomenta_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for momenta
// This implementation uses an AOSOA[npagM][npar][np4][neppM] where nevt=npagM*neppM
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessMomenta//_AOSOAv1
{
public:

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Indexed access (WITH an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  static
  __device__ inline
  fptype& ieventAccessIp4Ipar( fptype* buffer,
                               const int ievt,
                               const int ip4,
                               const int ipar )
  {
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    //printf( "%2d %2d %8d %8.3f\n", ipar, 0, ievt, buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
    return buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  // (Const memory access)
  static
  __device__ inline
  const fptype& ieventConstAccessIp4Ipar( fptype* buffer,
                                          const int ievt,
                                          const int ip4,
                                          const int ipar )
  {
    return ieventAccessIp4Ipar( const_cast<fptype*>( buffer ), ievt, ip4, ipar );
  }

};

//----------------------------------------------------------------------------

// A class describing host kernel access to memory buffers for momenta
class HostAccessMomenta : public MemoryAccessMomenta
{
public:

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Kernel access (WITHOUT an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  static
  __device__ inline
  fptype& kernelAccessIp4Ipar( fptype* buffer,
                               const int ip4,
                               const int ipar )
  {
    return ieventAccessIp4Ipar( buffer, 0, ip4, ipar );
  }

  // (Const memory access)
  static
  __device__ inline
  const fptype& kernelConstAccessIp4Ipar( fptype* buffer,
                                          const int ip4,
                                          const int ipar )
  {
    return kernelAccessIp4Ipar( const_cast<fptype*>( buffer ), ip4, ipar );
  }

};

//----------------------------------------------------------------------------

#ifdef __CUDACC__
// A class describing device kernel access to memory buffers for momenta
class DeviceAccessMomenta : public MemoryAccessMomenta
{
public:

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Kernel access (WITHOUT an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  static
  __device__ inline
  fptype& kernelAccessIp4Ipar( fptype* buffer,
                               const int ip4,
                               const int ipar )
  {
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "kernelAccessIp4Ipar: ievt=%d threadId=%d\n", ievt, threadIdx.x );
    return ieventAccessIp4Ipar( buffer, ievt, ip4, ipar ); // NB fptype and fptype_sv coincide for CUDA
  }

  // (Const memory access)
  static
  __device__ inline
  const fptype& kernelConstAccessIp4Ipar( fptype* buffer,
                                          const int ip4,
                                          const int ipar )
  {
    return kernelAccessIp4Ipar( const_cast<fptype*>( buffer ), ip4, ipar );
  }

};
#endif

//----------------------------------------------------------------------------

#endif // MemoryAccessMomenta_H
