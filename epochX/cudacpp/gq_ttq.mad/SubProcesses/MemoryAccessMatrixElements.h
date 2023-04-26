// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Jan 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2022-2023) for the MG5aMC CUDACPP plugin.

#ifndef MemoryAccessMatrixElements_H
#define MemoryAccessMatrixElements_H 1

#include "mgOnGpuConfig.h"

#include "MemoryAccessHelpers.h"
#include "MemoryAccessVectors.h"
#include "MemoryBuffers.h" // for HostBufferMatrixElements::isaligned

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for matrix elements
// This implementation uses a plain ARRAY[nevt]
// [If many implementations are used, a suffix _ARRAYv1 should be appended to the class name]
class MemoryAccessMatrixElementsBase //_ARRAYv1
{
private:

  friend class MemoryAccessHelper<MemoryAccessMatrixElementsBase>;
  friend class KernelAccessHelper<MemoryAccessMatrixElementsBase, true>;
  friend class KernelAccessHelper<MemoryAccessMatrixElementsBase, false>;

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
    return &( buffer[ievt] ); // ARRAY[nevt]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
  // [NB: expand variadic template "Ts... args" to empty and rename "Field" as empty]
  static __host__ __device__ inline fptype&
  decodeRecord( fptype* buffer )
  {
    constexpr int ievt = 0;
    return buffer[ievt]; // ARRAY[nevt]
  }
};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
// Its methods use the MemoryAccessHelper templates - note the use of the template keyword in template function instantiations
class MemoryAccessMatrixElements : public MemoryAccessMatrixElementsBase
{
public:

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (non-const) ===> fptype* ieventAccessRecord( fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::ieventAccessRecord;

  // Locate an event record (output) in a memory buffer (input) from the given event number (input)
  // [Signature (const) ===> const fptype* ieventAccessRecordConst( const fptype* buffer, const int ievt ) <===]
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::ieventAccessRecordConst;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer ) <===]
  static constexpr auto decodeRecord = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::decodeRecord;

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // [Signature (const) ===> const fptype& decodeRecordConst( const fptype* buffer ) <===]
  static constexpr auto decodeRecordConst =
    MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template decodeRecordConst<>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& ieventAccess( fptype* buffer, const ievt ) <===]
  static constexpr auto ieventAccess =
    MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template ieventAccessField<>;

  // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
  // [Signature (const) ===> const fptype& ieventAccessConst( const fptype* buffer, const ievt ) <===]
  static constexpr auto ieventAccessConst =
    MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template ieventAccessFieldConst<>;
};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
// Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
template<bool onDevice>
class KernelAccessMatrixElements
{
public:

  // Expose selected functions from MemoryAccessMatrixElements
  static constexpr auto ieventAccessRecord = MemoryAccessMatrixElements::ieventAccessRecord;

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const, SCALAR) ===> fptype& kernelAccess_s( fptype* buffer ) <===]
  static constexpr auto kernelAccess_s =
    KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessField<>; // requires cuda 11.4

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal)
  // [Signature (non const, SCALAR OR VECTOR) ===> fptype_sv& kernelAccess( const fptype* buffer ) <===]
  static __host__ __device__ inline fptype_sv&
  kernelAccess( fptype* buffer )
  {
    fptype& out = kernelAccess_s( buffer );
#ifndef MGONGPU_CPPSIMD
    return out;
#else
    // NB: derived from MemoryAccessMomenta, restricting the implementation to contiguous aligned arrays (#435)
    static_assert( mg5amcCpu::HostBufferMatrixElements::isaligned() ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
    //assert( (size_t)( buffer ) % mgOnGpu::cppAlign == 0 ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
    return mg5amcCpu::fptypevFromAlignedArray( out ); // SIMD bulk load of neppV, use reinterpret_cast
#endif
  }

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const) ===> const fptype& kernelAccessConst( const fptype* buffer ) <===]
  static constexpr auto kernelAccessConst =
    KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessFieldConst<>; // requires cuda 11.4
};

//----------------------------------------------------------------------------

typedef KernelAccessMatrixElements<false> HostAccessMatrixElements;
typedef KernelAccessMatrixElements<true> DeviceAccessMatrixElements;

//----------------------------------------------------------------------------

#endif // MemoryAccessMatrixElements_H
