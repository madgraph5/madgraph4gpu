#ifndef MemoryAccessMomenta_H
#define MemoryAccessMomenta_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

#include "MemoryAccessBase.h"

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
  // (Non-const memory access to event record from ievent)
  static
  __host__ __device__ inline
  fptype* ieventAccessRecord( fptype* buffer,
                              const int ievt )
  {
    constexpr int ip4 = 0;
    constexpr int ipar = 0;
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    return &( buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] ); // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // (Non-const memory access to field in an event record)
  static
  __host__ __device__ inline
  fptype& decodeRecord( fptype* buffer,
                        const int ip4,
                        const int ipar )
  {
    constexpr int ipagM = 0;
    constexpr int ieppM = 0;
    return buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  static constexpr auto decodeRecordIp4Ipar = decodeRecord;

  // *** BOILERPLATE STARTS ***

  // (Const memory access to event record from ievent)
  static
  __host__ __device__ inline
  fptype* ieventAccessConstRecord( const fptype* buffer,
                                   const int ievt )
  {
    return ieventAccessRecord( const_cast<fptype*>( buffer ), ievt );
  }

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

  // *** BOILERPLATE ENDS ***

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccessIp4Ipar = MemoryAccessBase::ieventAccessField<MemoryAccessMomenta>;

  // (Const memory access to field from ievent)
  static constexpr auto ieventConstAccessIp4Ipar = ieventAccessConstField;

};

//----------------------------------------------------------------------------

// A class describing kernel access to memory buffers for momenta on a CPU host or on a GPU device
template<bool onDevice>
class KernelAccessMomenta : public MemoryAccessMomenta
{
public:

  // *** BOILERPLATE STARTS ***

  // Locate an event record (output) in a memory buffer (input) from an implicit event-indexing mechanism in the kernel
  // (Non-const memory access to event record from kernel)
  static
  __host__ __device__ inline
  fptype* kernelAccessRecord( fptype* buffer )
  {
    //if constexpr ( !onDevice ) // FIXME! enable this when we move to nvcc supporting c++17
    if ( !onDevice )
    {
      return ieventAccessRecord( buffer, 0 );
    }
    else
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "kernelAccessRecord: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      return ieventAccessRecord( buffer, ievt ); // NB fptype and fptype_sv coincide for CUDA
#else
      throw std::runtime_error( "kernelAccessRecord on device is only implemented in CUDA" );
#endif
    }
  }

  // (Const memory access to event record from kernel)
  static
  __host__ __device__ inline
  fptype* kernelAccessConstRecord( const fptype* buffer )
  {
    return kernelAccessRecord( const_cast<fptype*>( buffer ) );
  }

  // (Non-const memory access to field from kernel)
  static
  __host__ __device__ inline
  fptype& kernelAccessField( fptype* buffer,
                             const int ip4,
                             const int ipar )
  {
    // NB all KernelLauncher classes assume that memory access can be decomposed in this way
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    return decodeRecord( kernelAccessRecord( buffer ), ip4, ipar );
  }

  // (Const memory access to field from ievent)
  static
  __host__ __device__ inline
  const fptype& kernelAccessConstField( const fptype* buffer,
                                        const int ip4,
                                        const int ipar )
  {
    return kernelAccessField( const_cast<fptype*>( buffer ), ip4, ipar );
  }

  // *** BOILERPLATE ENDS ***

  // (Non-const memory access to field from ievent)
  static constexpr auto kernelAccessIp4Ipar = kernelAccessField;

  // (Const memory access to field from ievent)
  static constexpr auto kernelConstAccessIp4Ipar = kernelAccessConstField;

};

//----------------------------------------------------------------------------

typedef KernelAccessMomenta<false> HostAccessMomenta;
typedef KernelAccessMomenta<true> DeviceAccessMomenta;

//----------------------------------------------------------------------------

#endif // MemoryAccessMomenta_H
