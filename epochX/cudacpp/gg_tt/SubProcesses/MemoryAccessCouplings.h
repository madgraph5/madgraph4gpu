#ifndef MemoryAccessCouplings_H
#define MemoryAccessCouplings_H 1

#include "mgOnGpuConfig.h"

#include "mgOnGpuCxtypes.h"

#include "MemoryAccessHelpers.h"
#include "MemoryAccessMomenta.h" // for MemoryAccessMomentaBase::neppM
#include "MemoryBuffers.h" // for HostBufferCouplings::isaligned

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for couplings
// This implementation uses an AOSOA[npagW][nx2][neppW] where nevt=npagW*neppW
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessCouplingsBase //_AOSOAv1
{
public:

  // Number of Events Per Page in the coupling AOSOA memory buffer layout
  static constexpr int neppC = MemoryAccessMomentaBase::neppM; // use the same AOSOA striding as for momenta

  // SANITY CHECK: check that neppC is a power of two
  static_assert( ispoweroftwo( neppC ), "neppC is not a power of 2" );

private:

  friend class MemoryAccessHelper<MemoryAccessCouplingsBase>;
  friend class KernelAccessHelper<MemoryAccessCouplingsBase, true>;
  friend class KernelAccessHelper<MemoryAccessCouplingsBase, false>;

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
    constexpr int ix2 = 0;
    const int ipagC = ievt / neppC;                                // #event "C-page"
    const int ieppC = ievt % neppC;                                // #event in the current event C-page
    return &( buffer[ipagC * nx2 * neppC + ix2 * neppC + ieppC] ); // AOSOA[ipagC][ix2][ieppC]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
  // [NB: expand variadic template "Ts... args" to "const int ix2" and rename "Field" as "Ix2"]
  static __host__ __device__ inline fptype&
  decodeRecord( fptype* buffer,
                const int ix2 )
  {
    constexpr int ipagC = 0;
    constexpr int ieppC = 0;
    return buffer[ipagC * nx2 * neppC + ix2 * neppC + ieppC]; // AOSOA[ipagC][ix2][ieppC]
  }
};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
// Its methods use the MemoryAccessHelper templates - note the use of the template keyword in template function instantiations
class MemoryAccessCouplings : public MemoryAccessCouplingsBase
{
public:

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessCouplingsBase>::ieventAccessRecord;

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (const) ===> const fptype* ieventAccessRecordConst( const fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessCouplingsBase>::ieventAccessRecordConst;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, const int ix2 ) <===]
  static constexpr auto decodeRecordIx2 = MemoryAccessHelper<MemoryAccessCouplingsBase>::decodeRecord;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (const) ===> const fptype& decodeRecordConst( const fptype* buffer, const int ix2 ) <===]
  static constexpr auto decodeRecordIx2Const =
    MemoryAccessHelper<MemoryAccessCouplingsBase>::template decodeRecordConst<int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& ieventAccessIx2( fptype* buffer, const ievt, const int ix2 ) <===]
  static constexpr auto ieventAccessIx2 =
    MemoryAccessHelper<MemoryAccessCouplingsBase>::template ieventAccessField<int>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (const) ===> const fptype& ieventAccessIx2Const( const fptype* buffer, const ievt, const int ix2 ) <===]
  static constexpr auto ieventAccessIx2Const =
    MemoryAccessHelper<MemoryAccessCouplingsBase>::template ieventAccessFieldConst<int>;
};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
// Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
template<bool onDevice>
class KernelAccessCouplings
{
public:

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const, SCALAR) ===> fptype& kernelAccessIx2( fptype* buffer, const int ix2 ) <===]
  static constexpr auto kernelAccessIx2_s =
    KernelAccessHelper<MemoryAccessCouplingsBase, onDevice>::template kernelAccessField<int>;

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const, SCALAR) ===> const fptype& kernelAccessIx2Const( const fptype* buffer, const int ix2 ) <===]
  static constexpr auto kernelAccessIx2Const_s =
    KernelAccessHelper<MemoryAccessCouplingsBase, onDevice>::template kernelAccessFieldConst<int>;

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non const, SCALAR OR VECTOR) ===> fptype_sv& kernelAccessIx2( fptype* buffer, const int ix2 ) <===]
  static __host__ __device__ inline fptype_sv&
  kernelAccessIx2( fptype* buffer,
                   const int ix2 )
  {
    fptype& out = kernelAccessIx2_s( buffer, ix2 );
#ifndef MGONGPU_CPPSIMD
    return out;
#else
    // NB: derived from MemoryAccessMomenta, restricting the implementation to contiguous aligned arrays
    constexpr int neppC = MemoryAccessCouplingsBase::neppC;
    static_assert( neppC >= neppV ); // ASSUME CONTIGUOUS ARRAYS
    static_assert( neppC % neppV == 0 ); // ASSUME CONTIGUOUS ARRAYS
    static_assert( mg5amcCpu::HostBufferCouplings::isaligned() ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)    
    //assert( (size_t)( buffer ) % mgOnGpu::cppAlign == 0 ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
    return mg5amcCpu::fptypevFromAlignedArray( out ); // SIMD bulk load of neppV, use reinterpret_cast
#endif
  }

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const, SCALAR OR VECTOR) ===> const fptype_sv& kernelAccessIx2Const( const fptype* buffer, const int ix2 ) <===]
  static __host__ __device__ inline const fptype_sv&
  kernelAccessIx2Const( const fptype* buffer,
                        const int ix2 )
  {
    return kernelAccessIx2Const( const_cast<fptype*>( buffer ), ix2 );
  }

  /*
  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const, SCALAR OR VECTOR) ===> const fptype_sv& kernelAccessIx2Const( const fptype* buffer, const int ix2 ) <===]
  static __host__ __device__ inline const fptype_sv&
  kernelAccessIx2Const( const fptype* buffer,
                        const int ix2 )
  {
    const fptype& out = kernelAccessIx2Const_s( buffer, ix2 );
#ifndef MGONGPU_CPPSIMD
    return out;
#else
    // NB: derived from MemoryAccessMomenta, restricting the implementation to contiguous aligned arrays
    constexpr int neppC = MemoryAccessCouplingsBase::neppC;
    static_assert( neppC >= neppV ); // ASSUME CONTIGUOUS ARRAYS
    static_assert( neppC % neppV == 0 ); // ASSUME CONTIGUOUS ARRAYS
    static_assert( mg5amcCpu::HostBufferCouplings::isaligned() ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)    
    //assert( (size_t)( buffer ) % mgOnGpu::cppAlign == 0 ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
    return mg5amcCpu::fptypevFromAlignedArray( out ); // SIMD bulk load of neppV, use reinterpret_cast
#endif
  }
  */
  
  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non const, SCALAR OR VECTOR) ===> cxtype_sv_ref kernelAccessConst( fptype* buffer ) <===]
  static __host__ __device__ inline cxtype_sv_ref
  kernelAccess( fptype* buffer )
  {
    return cxtype_sv_ref( kernelAccessIx2( buffer, 0 ),
                          kernelAccessIx2( buffer, 1 ) );
  }

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const, SCALAR OR VECTOR) ===> cxtype_sv kernelAccessConst( const fptype* buffer ) <===]
  static __host__ __device__ inline cxtype_sv
  kernelAccessConst( const fptype* buffer )
  {
    return cxtype_sv( kernelAccessIx2Const( buffer, 0 ),
                      kernelAccessIx2Const( buffer, 1 ) );
  }
};
  
//----------------------------------------------------------------------------

typedef KernelAccessCouplings<false> HostAccessCouplings;
typedef KernelAccessCouplings<true> DeviceAccessCouplings;

//----------------------------------------------------------------------------

#endif // MemoryAccessCouplings_H
