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
class MemoryAccessMomentaBase//_AOSOAv1
{
public:

  static constexpr int np4 = mgOnGpu::np4;
  static constexpr int npar = mgOnGpu::npar;
  static constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time

  //--------------------------------------------------------------------------

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

  //--------------------------------------------------------------------------

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

};

//----------------------------------------------------------------------------

//template MemoryAccessBase<MemoryAccessMomentaBase>::ieventAccessFIELD<const int, const int>;
//template fptype& MemoryAccessBase<MemoryAccessMomentaBase>::ieventAccessFIELD<int, int>( fptype*, int, int, int );

// A class providing access to memory buffers for a given event, based on explicit event numbers
class MemoryAccessMomenta : public MemoryAccessMomentaBase
{
public:

  // (Non-const memory access to field in an event record)
  //static constexpr auto decodeRecordIp4Ipar = MemoryAccessBase<MemoryAccessMomentaBase>::decodeRecord;

  // (Non-const memory access to field from ievent)
  //static constexpr auto ieventAccessIp4Ipar = MemoryAccessBase<MemoryAccessMomentaBase>::ieventAccessField;
  static constexpr auto ieventAccessIp4Ipar = MemoryAccessBase<MemoryAccessMomentaBase>::ieventAccessFIELD<int, int>;
  //typedef fptype& (*functype) (fptype*, const int, const int, const int);
  //static constexpr functype ieventAccessIp4Ipar = MemoryAccessBase<MemoryAccessMomentaBase>::ieventAccessFIELD<const int, const int>;

  // (Const memory access to field from ievent)
  //static constexpr auto ieventConstAccessIp4Ipar = MemoryAccessBase<MemoryAccessMomentaBase>::ieventAccessConstField;

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
template<bool onDevice>
class KernelAccessMomenta
{
public:

  // (Non-const memory access to field from ievent)
  static constexpr auto kernelAccessIp4Ipar = KernelAccessBase<MemoryAccessMomentaBase, onDevice>::kernelAccessField;

  // (Const memory access to field from ievent)
  //static constexpr auto kernelConstAccessIp4Ipar = KernelAccessBase<MemoryAccessMomentaBase, onDevice>::kernelAccessConstField;

};

//----------------------------------------------------------------------------

typedef KernelAccessMomenta<false> HostAccessMomenta;
typedef KernelAccessMomenta<true> DeviceAccessMomenta;

//----------------------------------------------------------------------------

#endif // MemoryAccessMomenta_H
