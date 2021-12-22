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

  static constexpr int np4 = mgOnGpu::np4;
  static constexpr int npar = mgOnGpu::npar;
  static constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time

  // Locate an event record (output) in a memory buffer (input) from an explicit event number (input)
  // [Rename as getRecordIevent?]
  // (Non-const memory access)
  static
  __host__ __device__ inline
  fptype* ieventAccess( fptype* buffer,
                        const int ievt )
  {
    constexpr int ip4 = 0;
    constexpr int ipar = 0;
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    return &( buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] ); // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Rename as decodeRecordIp4Ipar?]
  // (Non-const memory access)
  static
  __host__ __device__ inline
  fptype& decodeRecordIp4Ipar( fptype* buffer,
                               const int ip4,
                               const int ipar )
  {
    constexpr int ipagM = 0;
    constexpr int ieppM = 0;
    return buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  // *** EVERYTHING BELOW IS BOILERPLATE ***

  // =========================================================================
  // *** Pattern: ieventAccessInd1..IndN( buffer, ievt [, ind1... indN] )  ***
  // =========================================================================

  // Indexed access (WITH an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  static
  __host__ __device__ inline
  fptype& ieventAccessIp4Ipar( fptype* buffer,
                               const int ievt,
                               const int ip4,
                               const int ipar )
  {
    // NB all KernelLauncher classes assume that memory access can be decomposed in this way
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return decodeRecordIp4Ipar( ieventAccess( buffer, ievt ), ip4, ipar );
  }

  // (Const memory access)
  static
  __host__ __device__ inline
  const fptype& ieventConstAccessIp4Ipar( fptype* buffer,
                                          const int ievt,
                                          const int ip4,
                                          const int ipar )
  {
    return ieventAccessIp4Ipar( const_cast<fptype*>( buffer ), ievt, ip4, ipar );
  }

};

//----------------------------------------------------------------------------

// A class describing kernel access to memory buffers for momenta on a CPU host or on a GPU device
template<bool onDevice>
class KernelAccessMomenta : public MemoryAccessMomenta
{
public:

  // =========================================================================
  // *** Pattern: kernelAccessInd1..IndN( buffer [, ind1... indN] )        ***
  // =========================================================================

  // Kernel access (WITHOUT an explicit event number) to the memory buffer for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: a specific 4-momenta component for a specific particle in one event, given its event number
  // (Non-const memory access)
  static
  __host__ __device__ inline
  fptype& kernelAccessIp4Ipar( fptype* buffer,
                               const int ip4,
                               const int ipar )
  {
    //if constexpr ( !onDevice ) // FIXME! enable this when we move to nvcc supporting c++17
    if ( !onDevice )
    {
      return ieventAccessIp4Ipar( buffer, 0, ip4, ipar );
    }
    else
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "kernelAccessIp4Ipar: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      return ieventAccessIp4Ipar( buffer, ievt, ip4, ipar ); // NB fptype and fptype_sv coincide for CUDA
#else
      throw std::runtime_error( "KernelAccessMomenta on device is only implemented in CUDA" );
#endif
    }
  }

  // (Const memory access)
  static
  __host__ __device__ inline
  const fptype& kernelConstAccessIp4Ipar( fptype* buffer,
                                          const int ip4,
                                          const int ipar )
  {
    return kernelAccessIp4Ipar( const_cast<fptype*>( buffer ), ip4, ipar );
  }

};

//----------------------------------------------------------------------------

typedef KernelAccessMomenta<false> HostAccessMomenta;
typedef KernelAccessMomenta<true> DeviceAccessMomenta;

//----------------------------------------------------------------------------

#endif // MemoryAccessMomenta_H
