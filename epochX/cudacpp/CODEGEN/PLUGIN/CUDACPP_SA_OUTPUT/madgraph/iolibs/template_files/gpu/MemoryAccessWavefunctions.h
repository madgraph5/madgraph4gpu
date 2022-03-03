#ifndef MemoryAccessWavefunctions_H
#define MemoryAccessWavefunctions_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "MemoryAccessHelpers.h"

#define MGONGPU_TRIVIAL_WAVEFUNCTIONS 1

//----------------------------------------------------------------------------

#ifndef MGONGPU_TRIVIAL_WAVEFUNCTIONS

// A class describing the internal layout of memory buffers for wavefunctions
// This implementation uses an AOSOA[npagW][nw6][nx2][neppW] where nevt=npagW*neppW
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessWavefunctionsBase //_AOSOAv1
{
public:

  // Number of Events Per Page in the wavefunction AOSOA memory buffer layout
  static constexpr int neppW = 1; // AOS (just a test...)

private:

  friend class MemoryAccessHelper<MemoryAccessWavefunctionsBase>;
  friend class KernelAccessHelper<MemoryAccessWavefunctionsBase, true>;
  friend class KernelAccessHelper<MemoryAccessWavefunctionsBase, false>;

  // The number of components of a (fermion or vector) wavefunction
  static constexpr int nw6 = mgOnGpu::nw6;

  // The number of floating point components of a complex number
  static constexpr int nx2 = mgOnGpu::nx2;

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
    const int ipagW = ievt / neppW; // #event "W-page"
    const int ieppW = ievt % neppW; // #event in the current event W-page
    constexpr int iw6 = 0;
    constexpr int ix2 = 0;
    return &( buffer[ipagW * nw6 * nx2 * neppW + iw6 * nx2 * neppW + ix2 * neppW + ieppW] ); // AOSOA[ipagW][iw6][ix2][ieppW]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
  // [NB: expand variadic template "Ts... args" to "const int ix2, const int iw6" and rename "Field" as "Ix2Iw6"]
  static __host__ __device__ inline fptype&
  decodeRecord( fptype* buffer,
                const int ix2,
                const int iw6 )
  {
    constexpr int ipagW = 0;
    constexpr int ieppW = 0;
    return buffer[ipagW * nw6 * nx2 * neppW + iw6 * nx2 * neppW + ix2 * neppW + ieppW]; // AOSOA[ipagW][iw6][ix2][ieppW]
  }
};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
// Its methods use the MemoryAccessHelper templates - note the use of the template keyword in template function instantiations
class MemoryAccessWavefunctions : public MemoryAccessWavefunctionsBase
{
public:

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessWavefunctionsBase>::ieventAccessRecord;

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (const) ===> const fptype* ieventAccessRecordConst( const fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessWavefunctionsBase>::ieventAccessRecordConst;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, const int ipar, const int iw6 ) <===]
  static constexpr auto decodeRecordIx2Iw6 = MemoryAccessHelper<MemoryAccessWavefunctionsBase>::decodeRecord;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (const) ===> const fptype& decodeRecordConst( const fptype* buffer, const int ipar, const int iw6 ) <===]
  static constexpr auto decodeRecordIx2Iw6Const =
    MemoryAccessHelper<MemoryAccessWavefunctionsBase>::template decodeRecordConst<int, int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& ieventAccessIx2Iw6( fptype* buffer, const ievt, const int ipar, const int iw6 ) <===]
  static constexpr auto ieventAccessIx2Iw6 =
    MemoryAccessHelper<MemoryAccessWavefunctionsBase>::template ieventAccessField<int, int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (const) ===> const fptype& ieventAccessIx2Iw6Const( const fptype* buffer, const ievt, const int ipar, const int iw6 ) <===]
  static constexpr auto ieventAccessIx2Iw6Const =
    MemoryAccessHelper<MemoryAccessWavefunctionsBase>::template ieventAccessFieldConst<int, int>;
};

#endif // #ifndef MGONGPU_TRIVIAL_WAVEFUNCTIONS

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
// Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
template<bool onDevice>
class KernelAccessWavefunctions
{
public:

#ifndef MGONGPU_TRIVIAL_WAVEFUNCTIONS

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& kernelAccessIx2Iw6( fptype* buffer, const int ipar, const int iw6 ) <===]
  static constexpr auto kernelAccessIx2Iw6 =
    KernelAccessHelper<MemoryAccessWavefunctionsBase, onDevice>::template kernelAccessField<int, int>;

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const) ===> const fptype& kernelAccessIx2Iw6Const( const fptype* buffer, const int ipar, const int iw6 ) <===]
  static constexpr auto kernelAccessIx2Iw6Const =
    KernelAccessHelper<MemoryAccessWavefunctionsBase, onDevice>::template kernelAccessFieldConst<int, int>;

#else

  static __host__ __device__ inline cxtype_sv*
  kernelAccess( fptype* buffer )
  {
    return reinterpret_cast<cxtype_sv*>( buffer );
  }

  static __host__ __device__ inline const cxtype_sv*
  kernelAccessConst( const fptype* buffer )
  {
    return reinterpret_cast<const cxtype_sv*>( buffer );
  }

#endif // #ifndef MGONGPU_TRIVIAL_WAVEFUNCTIONS
};

//----------------------------------------------------------------------------

typedef KernelAccessWavefunctions<false> HostAccessWavefunctions;
typedef KernelAccessWavefunctions<true> DeviceAccessWavefunctions;

//----------------------------------------------------------------------------

#endif // MemoryAccessWavefunctions_H
