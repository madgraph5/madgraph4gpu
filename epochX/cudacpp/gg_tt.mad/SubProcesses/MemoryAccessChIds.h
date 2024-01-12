// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: S. Roiser (Dec 2023) for the MG5aMC CUDACPP plugin.

#ifndef MemoryAccessChIds_H
#define MemoryAccessChIds_H 1

#include "mgOnGpuConfig.h"

#include "MemoryAccessHelpers.h"
#include "MemoryAccessVectors.h"
#include "MemoryBuffers.h" // for HostBufferMatrixElements::isaligned

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //----------------------------------------------------------------------------

  // A class describing the internal layout of memory buffers for Channel Ids
  // This implementation uses a plain ARRAY[nevt]
  // [If many implementations are used, a suffix _ARRAYv1 should be appended to the class name]
  class MemoryAccessChIdsBase //_ARRAYv1
  {
  private:

    friend class MemoryAccessHelper<MemoryAccessChIdsBase>;
    friend class KernelAccessHelper<MemoryAccessChIdsBase, true>;
    friend class KernelAccessHelper<MemoryAccessChIdsBase, false>;

    //--------------------------------------------------------------------------
    // NB all KernelLaunchers assume that memory access can be decomposed as "accessField = decodeRecord( accessRecord )"
    // (in other words: first locate the event record for a given event, then locate an element in that record)
    //--------------------------------------------------------------------------

    // Locate an event record (output) in a memory buffer (input) from the given event number (input)
    // [Signature (non-const) ===> unsigned int* ieventAccessRecord( unsigned int* buffer, const int ievt ) <===]
    static __host__ __device__ inline unsigned int*
    ieventAccessRecord( unsigned int* buffer,
                        const int ievt )
    {
      return &( buffer[ievt] ); // ARRAY[nevt]
    }

    //--------------------------------------------------------------------------

    // Locate a field (output) of an event record (input) from the given field indexes (input)
    // [Signature (non-const) ===> unsigned int& decodeRecord( unsigned int* buffer, Ts... args ) <===]
    // [NB: expand variadic template "Ts... args" to empty and rename "Field" as empty]
    static __host__ __device__ inline unsigned int&
    decodeRecord( unsigned int* buffer )
    {
      constexpr int ievt = 0;
      return buffer[ievt]; // ARRAY[nevt]
    }
  };

  //----------------------------------------------------------------------------

  // A class providing access to memory buffers for a given event, based on explicit event numbers
  // Its methods use the MemoryAccessHelper templates - note the use of the template keyword in template function instantiations
  class MemoryAccessChIds : public MemoryAccessChIdsBase
  {
  public:

    // Locate an event record (output) in a memory buffer (input) from the given event number (input)
    // [Signature (non-const) ===> unsigned int* ieventAccessRecord( unsigned int* buffer, const int ievt ) <===]
    static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessChIdsBase>::ieventAccessRecord;

    // Locate an event record (output) in a memory buffer (input) from the given event number (input)
    // [Signature (const) ===> const unsigned int* ieventAccessRecordConst( const unsigned int* buffer, const int ievt ) <===]
    static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessChIdsBase>::ieventAccessRecordConst;

    // Locate a field (output) of an event record (input) from the given field indexes (input)
    // [Signature (non-const) ===> unsigned int& decodeRecord( unsigned int* buffer ) <===]
    static constexpr auto decodeRecord = MemoryAccessHelper<MemoryAccessChIdsBase>::decodeRecord;

    // Locate a field (output) of an event record (input) from the given field indexes (input)
    // [Signature (const) ===> const unsigned int& decodeRecordConst( const unsigned int* buffer ) <===]
    static constexpr auto decodeRecordConst =
      MemoryAccessHelper<MemoryAccessChIdsBase>::template decodeRecordConst<>;

    // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
    // [Signature (non-const) ===> unsigned int& ieventAccess( unsigned int* buffer, const ievt ) <===]
    static constexpr auto ieventAccess =
      MemoryAccessHelper<MemoryAccessChIdsBase>::template ieventAccessField<>;

    // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
    // [Signature (const) ===> const unsigned int& ieventAccessConst( const unsigned int* buffer, const ievt ) <===]
    static constexpr auto ieventAccessConst =
      MemoryAccessHelper<MemoryAccessChIdsBase>::template ieventAccessFieldConst<>;
  };

  //----------------------------------------------------------------------------

  // A class providing access to memory buffers for a given event, based on implicit kernel rules
  // Its methods use the KernelAccessHelper template - note the use of the template keyword in template function instantiations
  template<bool onDevice>
  class KernelAccessChIds
  {
  public:

    // Expose selected functions from MemoryAccessChIds
    static constexpr auto ieventAccessRecord = MemoryAccessChIds::ieventAccessRecord;

    // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
    // [Signature (non-const, SCALAR) ===> unsigned int& kernelAccess( unsigned int* buffer ) <===]
    static constexpr auto kernelAccess_s =
      KernelAccessHelper<MemoryAccessChIdsBase, onDevice>::template kernelAccessField<>; // requires cuda 11.4

    // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal)
    // [Signature (non-const, SCALAR OR VECTOR) ===> fptype_sv& kernelAccess( fptype* buffer ) <===]
    static __host__ __device__ inline fptype_sv&
    kernelAccess( fptype* buffer )
    {
      fptype& out = kernelAccess_s( buffer );
#ifndef MGONGPU_CPPSIMD
      return out;
#else
      // NB: derived from MemoryAccessMomenta, restricting the implementation to contiguous aligned arrays (#435)
      static_assert( mg5amcCpu::HostBufferChIds::isaligned() ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
      //assert( (size_t)( buffer ) % mgOnGpu::cppAlign == 0 ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
      return mg5amcCpu::fptypevFromAlignedArray( out ); // SIMD bulk load of neppV, use reinterpret_cast
#endif
    }

    // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
    // [Signature (const, SCALAR) ===> const fptype& kernelAccessConst( const fptype* buffer ) <===]
    static constexpr auto kernelAccessConst_s =
      KernelAccessHelper<MemoryAccessChIdsBase, onDevice>::template kernelAccessFieldConst<>; // requires cuda 11.4

    // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal)
    // [Signature (const, SCALAR OR VECTOR) ===> const fptype_sv& kernelAccess( const fptype* buffer ) <===]
    static __host__ __device__ inline const fptype_sv&
    kernelAccessConst( const fptype* buffer )
    {
      const fptype& out = kernelAccessConst_s( buffer );
#ifndef MGONGPU_CPPSIMD
      return out;
#else
      // NB: derived from MemoryAccessMomenta, restricting the implementation to contiguous aligned arrays (#435)
      static_assert( mg5amcCpu::HostBufferChIds::isaligned() ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
      //assert( (size_t)( buffer ) % mgOnGpu::cppAlign == 0 ); // ASSUME ALIGNED ARRAYS (reinterpret_cast will segfault otherwise!)
      return mg5amcCpu::fptypevFromAlignedArray( out ); // SIMD bulk load of neppV, use reinterpret_cast
#endif
    }
  };

  //----------------------------------------------------------------------------

  typedef KernelAccessChIds<false> HostAccessChIds;
  typedef KernelAccessChIds<true> DeviceAccessChIds;

  //----------------------------------------------------------------------------

} // end namespace mg5amcGpu/mg5amcCpu

#endif // MemoryAccessChIds_H
