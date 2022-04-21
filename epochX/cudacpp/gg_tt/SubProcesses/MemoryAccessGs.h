#ifndef MemoryAccessGs_H
#define MemoryAccessGs_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "MemoryAccessHelpers.h"

#define MGONGPU_TRIVIAL_GS 1

//----------------------------------------------------------------------------

#ifndef MGONGPU_TRIVIAL_GS

// A class describing the internal layout of memory buffers for Gs
// This implementation uses an AOSOA[npagW][nx2][neppW] where nevt=npagW*neppW
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessGsBase //_AOSOAv1
{
public:

  // Number of Events Per Page in the G AOSOA memory buffer layout
  static constexpr int neppW = 1; // AOS (just a test...)

private:

  friend class MemoryAccessHelper<MemoryAccessGsBase>;
  friend class KernelAccessHelper<MemoryAccessGsBase, true>;
  friend class KernelAccessHelper<MemoryAccessGsBase, false>;

  // The number of floating point components of a complex number
  static constexpr int nx2 = 1; // mgOnGpu::nx2;

  //--------------------------------------------------------------------------
  // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
  // (in other words: first locate the event record for a given event, then locate an element in that record)
  //--------------------------------------------------------------------------

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static __host__ __device__ inline fptype*
  ieventAccessRecord( fptype* buffer,
                      const int ievt )
  {
    constexpr int ix2 = 0;
    const int ipagW = ievt / neppW;                                // #event "W-page"
    const int ieppW = ievt % neppW;                                // #event in the current event W-page
    return &( buffer[ipagW * nx2 * neppW + ix2 * neppW + ieppW] ); // AOSOA[ipagW][ix2][ieppW]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
  // [NB: expand variadic template "Ts... args" to "const int ix2" and rename "Field" as "Ix2"]
  static __host__ __device__ inline fptype&
  decodeRecord( fptype* buffer,
                const int ix2 )
  {
    constexpr int ipagW = 0;
    constexpr int ieppW = 0;
    return buffer[ipagW * nx2 * neppW + ix2 * neppW + ieppW]; // AOSOA[ipagW][ix2][ieppW]
  }
};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
// Its methods use the MemoryAccessHelper templates - note the use of the template keyword in template function instantiations
class MemoryAccessGs : public MemoryAccessGsBase
{
public:

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessGsBase>::ieventAccessRecord;

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (const) ===> const fptype* ieventAccessRecordConst( const fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessGsBase>::ieventAccessRecordConst;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, const int ix2 ) <===]
  static constexpr auto decodeRecordIx2 = MemoryAccessHelper<MemoryAccessGsBase>::decodeRecord;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (const) ===> const fptype& decodeRecordConst( const fptype* buffer, const int ix2 ) <===]
  static constexpr auto decodeRecordIx2Const =
    MemoryAccessHelper<MemoryAccessGsBase>::template decodeRecordConst<int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& ieventAccessIx2( fptype* buffer, const ievt, const int ix2 ) <===]
  static constexpr auto ieventAccessIx2 =
    MemoryAccessHelper<MemoryAccessGsBase>::template ieventAccessField<int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (const) ===> const fptype& ieventAccessIx2Const( const fptype* buffer, const ievt, const int ix2 ) <===]
  static constexpr auto ieventAccessIx2Const =
    MemoryAccessHelper<MemoryAccessGsBase>::template ieventAccessFieldConst<int>;
};

#endif // #ifndef MGONGPU_TRIVIAL_GS

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
// Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
template<bool onDevice>
class KernelAccessGs
{
public:

#ifndef MGONGPU_TRIVIAL_GS

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& kernelAccessIx2( fptype* buffer, const int ix2 ) <===]
  static constexpr auto kernelAccessIx2 =
    KernelAccessHelper<MemoryAccessGsBase, onDevice>::template kernelAccessField<int>;

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const) ===> const fptype& kernelAccessIx2Const( const fptype* buffer, const int ix2 ) <===]
  static constexpr auto kernelAccessIx2Const =
    KernelAccessHelper<MemoryAccessGsBase, onDevice>::template kernelAccessFieldConst<int>;

#else

  static __host__ __device__ inline fptype_sv*
  kernelAccess( fptype* buffer )
  {
    return reinterpret_cast<fptype_sv*>( buffer );
  }

  static __host__ __device__ inline const fptype_sv*
  kernelAccessConst( const fptype* buffer )
  {
    return reinterpret_cast<const fptype_sv*>( buffer );
  }

#endif // #ifndef MGONGPU_TRIVIAL_GS
};

//----------------------------------------------------------------------------

typedef KernelAccessGs<false> HostAccessGs;
typedef KernelAccessGs<true> DeviceAccessGs;

//----------------------------------------------------------------------------

#endif // MemoryAccessGs_H
