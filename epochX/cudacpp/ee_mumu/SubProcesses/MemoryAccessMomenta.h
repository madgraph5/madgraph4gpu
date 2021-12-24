#ifndef MemoryAccessMomenta_H
#define MemoryAccessMomenta_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

#include "MemoryAccessHelpers.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for momenta
// This implementation uses an AOSOA[npagM][npar][np4][neppM] where nevt=npagM*neppM
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessMomentaBase//_AOSOAv1
{
public:

  // Number of Events Per Page in the momenta AOSOA memory buffer layout
  static constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time

private:

  friend class MemoryAccessHelper<MemoryAccessMomentaBase>;
  friend class KernelAccessHelper<MemoryAccessMomentaBase, true>;
  friend class KernelAccessHelper<MemoryAccessMomentaBase, false>;
  
  // The number of components of a 4-momentum
  static constexpr int np4 = mgOnGpu::np4;

  // The number of particles in this physics process
  static constexpr int npar = mgOnGpu::npar;

  //--------------------------------------------------------------------------
  // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
  // (in other words: first locate the event record for a given event, then locate an element in that record)
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

// A class providing access to memory buffers for a given event, based on explicit event numbers
class MemoryAccessMomenta : public MemoryAccessMomentaBase
{
public:

  // (Non-const memory access to event record from ievent)
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessMomentaBase>::ieventAccessRecord;

  // (Const memory access to event record from ievent)
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessMomentaBase>::ieventAccessRecordConst;

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecordIp4Ipar = MemoryAccessHelper<MemoryAccessMomentaBase>::decodeRecord;

  // [NOTE THE USE OF THE TEMPLATE KEYWORD IN ALL OF THE FOLLOWING TEMPLATE FUNCTION INSTANTIATIONS]
  // (Const memory access to field in an event record)
  static constexpr auto decodeRecordIp4IparConst =
    MemoryAccessHelper<MemoryAccessMomentaBase>::template decodeRecordConst<int, int>;

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccessIp4Ipar =
    MemoryAccessHelper<MemoryAccessMomentaBase>::template ieventAccessField<int, int>;

#if 1
  // (Const memory access to field from ievent)
  static constexpr auto ieventAccessIp4IparConst =
    MemoryAccessHelper<MemoryAccessMomentaBase>::template ieventAccessFieldConst<int, int>;
#else
  // (Const memory access to field from ievent - DEBUG version with printouts)
  static
  __host__ __device__ inline
  const fptype& ieventAccessIp4IparConst( const fptype* buffer,
                                          const int ievt,
                                          const int ip4,
                                          const int ipar )
  {
    const fptype& out = MemoryAccessHelper<MemoryAccessMomentaBase>::template ieventAccessFieldConst<int, int>( buffer, ievt, ip4, ipar );
    printf( "ipar=%2d ip4=%2d ievt=%8d out=%8.3f\n", ipar, ip4, ievt, out );
    return out;
  }
#endif

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
template<bool onDevice>
class KernelAccessMomenta
{
public:

  // (Non-const memory access to field from kernel)
  static constexpr auto kernelAccessIp4Ipar =
    KernelAccessHelper<MemoryAccessMomentaBase, onDevice>::template kernelAccessField<int, int>;

#if 1
  // (Const memory access to field from kernel)
  static constexpr auto kernelAccessIp4IparConst =
    KernelAccessHelper<MemoryAccessMomentaBase, onDevice>::template kernelAccessFieldConst<int, int>;
#else
  // (Const memory access to field from kernel - DEBUG version with printouts)
  static
  __host__ __device__ inline
  const fptype& kernelAccessIp4IparConst( const fptype* buffer,
                                          const int ip4,
                                          const int ipar )
  {
    const fptype& out = KernelAccessHelper<MemoryAccessMomentaBase, onDevice>::template kernelAccessFieldConst<int, int>( buffer, ip4, ipar );
    printf( "ipar=%2d ip4=%2d ievt=  kernel out=%8.3f\n", ipar, ip4, out );
    return out;
  }
#endif

};

//----------------------------------------------------------------------------

typedef KernelAccessMomenta<false> HostAccessMomenta;
typedef KernelAccessMomenta<true> DeviceAccessMomenta;

//----------------------------------------------------------------------------

#endif // MemoryAccessMomenta_H
