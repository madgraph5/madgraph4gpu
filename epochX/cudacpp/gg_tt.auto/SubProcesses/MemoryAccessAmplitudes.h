#ifndef MemoryAccessAmplitudes_H
#define MemoryAccessAmplitudes_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuCxtypes.h"

#include "MemoryAccessHelpers.h"

#define MGONGPU_TRIVIAL_AMPLITUDES 1

//----------------------------------------------------------------------------

#ifndef MGONGPU_TRIVIAL_AMPLITUDES

// A class describing the internal layout of memory buffers for amplitudes
// This implementation uses an AOSOA[npagW][nx2][neppW] where nevt=npagW*neppW
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessAmplitudesBase//_AOSOAv1
{
public:

  // Number of Events Per Page in the amplitude AOSOA memory buffer layout
  static constexpr int neppW = 1; // AOS (just a test...)

 private:

  friend class MemoryAccessHelper<MemoryAccessAmplitudesBase>;
  friend class KernelAccessHelper<MemoryAccessAmplitudesBase, true>;
  friend class KernelAccessHelper<MemoryAccessAmplitudesBase, false>;
  
  // The number of floating point components of a complex number
  static constexpr int nx2 = mgOnGpu::nx2;

  //--------------------------------------------------------------------------
  // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
  // (in other words: first locate the event record for a given event, then locate an element in that record)
  //--------------------------------------------------------------------------

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static
  __host__ __device__ inline
  fptype* ieventAccessRecord( fptype* buffer,
                              const int ievt )
  {
    constexpr int ix2 = 0;
    const int ipagW = ievt/neppW; // #event "W-page"
    const int ieppW = ievt%neppW; // #event in the current event W-page
    return &( buffer[ipagW*nx2*neppW + ix2*neppW + ieppW] ); // AOSOA[ipagW][ix2][ieppW]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
  // [NB: expand variadic template "Ts... args" to "const int ix2" and rename "Field" as "Ix2"]
  static
  __host__ __device__ inline
  fptype& decodeRecord( fptype* buffer,
                        const int ix2 )
  {
    constexpr int ipagW = 0;
    constexpr int ieppW = 0;
    return buffer[ipagW*nx2*neppW + ix2*neppW + ieppW]; // AOSOA[ipagW][ix2][ieppW]
  }

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
// Its methods use the MemoryAccessHelper templates - note the use of the template keyword in template function instantiations
class MemoryAccessAmplitudes : public MemoryAccessAmplitudesBase
{
public:

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessAmplitudesBase>::ieventAccessRecord;

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (const) ===> const fptype* ieventAccessRecordConst( const fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessAmplitudesBase>::ieventAccessRecordConst;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, const int ix2 ) <===]
  static constexpr auto decodeRecordIx2 = MemoryAccessHelper<MemoryAccessAmplitudesBase>::decodeRecord;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (const) ===> const fptype& decodeRecordConst( const fptype* buffer, const int ix2 ) <===]
  static constexpr auto decodeRecordIx2Const =
    MemoryAccessHelper<MemoryAccessAmplitudesBase>::template decodeRecordConst<int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& ieventAccessIx2( fptype* buffer, const ievt, const int ix2 ) <===]
  static constexpr auto ieventAccessIx2 =
    MemoryAccessHelper<MemoryAccessAmplitudesBase>::template ieventAccessField<int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (const) ===> const fptype& ieventAccessIx2Const( const fptype* buffer, const ievt, const int ix2 ) <===]
  static constexpr auto ieventAccessIx2Const =
    MemoryAccessHelper<MemoryAccessAmplitudesBase>::template ieventAccessFieldConst<int>;

};

#endif // #ifndef MGONGPU_TRIVIAL_AMPLITUDES

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
// Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
template<bool onDevice>
class KernelAccessAmplitudes
{
public:

#ifndef MGONGPU_TRIVIAL_AMPLITUDES

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& kernelAccessIx2( fptype* buffer, const int ix2 ) <===]
  static constexpr auto kernelAccessIx2 =
    KernelAccessHelper<MemoryAccessAmplitudesBase, onDevice>::template kernelAccessField<int>;

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const) ===> const fptype& kernelAccessIx2Const( const fptype* buffer, const int ix2 ) <===]
  static constexpr auto kernelAccessIx2Const =
    KernelAccessHelper<MemoryAccessAmplitudesBase, onDevice>::template kernelAccessFieldConst<int>;

#else

  static
  __host__ __device__ inline
  cxtype_sv* kernelAccess( fptype* buffer )
  {
    return reinterpret_cast<cxtype_sv*>( buffer );
  }

  static
  __host__ __device__ inline
  const cxtype_sv* kernelAccessConst( const fptype* buffer )
  {
    return reinterpret_cast<const cxtype_sv*>( buffer );
  }

#endif // #ifndef MGONGPU_TRIVIAL_AMPLITUDES

};

//----------------------------------------------------------------------------

typedef KernelAccessAmplitudes<false> HostAccessAmplitudes;
typedef KernelAccessAmplitudes<true> DeviceAccessAmplitudes;

//----------------------------------------------------------------------------

#endif // MemoryAccessAmplitudes_H
