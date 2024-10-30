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

// NB: namespaces mg5amcGpu and mg5amcCpu includes types which are defined in different ways for CPU and GPU builds (see #318 and #725)
#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //----------------------------------------------------------------------------

  // A class describing the internal layout of memory buffers for matrix elements
  // This implementation uses an AOSOA[npagME][ncomb+1][neppME] where nevt=npagME*neppME
  // [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
  class MemoryAccessMatrixElementsBase //_AOSOAv1
  {
  public:

    // Number of Events Per Page in the matrix element AOSOA memory buffer layout
#ifdef MGONGPUCPP_GPUIMPL
    static constexpr int neppME = 32/sizeof(fptype); // (DEFAULT) 32-byte GPU cache line (256 bits): 4 (DOUBLE) or 8 (FLOAT)
#else
#ifdef MGONGPU_CPPSIMD
    static constexpr int neppME = MGONGPU_CPPSIMD; // (DEFAULT) neppME=neppV for optimal performance
#else
    static constexpr int neppME = 1; // (DEFAULT) neppM=neppV for optimal performance (NB: this is equivalent to AOS)
#endif
#endif

  private:

    friend class MemoryAccessHelper<MemoryAccessMatrixElementsBase>;
    friend class KernelAccessHelper<MemoryAccessMatrixElementsBase, true>;
    friend class KernelAccessHelper<MemoryAccessMatrixElementsBase, false>;

    // The number of (good and bad) helicity combinations
    // Note that a typical loop over helicities uses ihel over ncomb values: "for( int ihel = 0; ihel < ncomb; ihel++ )"
    static constexpr int ncomb = CPPProcess::ncomb;

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
      const int ipagME = ievt / neppME; // #event "ME-page"
      const int ieppME = ievt % neppME; // #event in the current event ME-page
      constexpr int ihel = 0;
      return &( buffer[ipagME * ( ncomb + 1 ) * neppME + ihel * neppME + ieppME] ); // AOSOA[ipagME][ihel][ieppME]
    }

    //--------------------------------------------------------------------------

    // Locate a field (output) of an event record (input) from the given field indexes (input)
    // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, Ts... args ) <===]
    // [NB: expand variadic template "Ts... args" to "const int ihel" and rename "Field" as "Ihel"]
    static __host__ __device__ inline fptype&
    decodeRecord( fptype* buffer,
                  const int ihel )
    {
      constexpr int ipagME = 0;
      constexpr int ieppME = 0;
      return buffer[ipagME * ( ncomb + 1 ) * neppME + ihel * neppME + ieppME]; // AOSOA[ipagME][ihel][ieppME]
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
    // [Signature (non-const) ===> fptype& decodeRecord( fptype* buffer, const int ihel ) <===]
    static constexpr auto decodeRecordIhel = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::decodeRecord;

    // Locate a field (output) of an event record (input) from the given field indexes (input)
    // [Signature (const) ===> const fptype& decodeRecordConst( const fptype* buffer, const int ihel ) <===]
    static constexpr auto decodeRecordIhelConst =
      MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template decodeRecordConst<int>;

    // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
    // [Signature (non-const) ===> fptype& ieventAccess( fptype* buffer, const ievt, const ihel ) <===]
    static constexpr auto ieventAccessIhel =
      MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template ieventAccessField<int>;

    // Locate a field (output) in a memory buffer (input) from the given event number (input) and the given field indexes (input)
    // [Signature (const) ===> const fptype& ieventAccessConst( const fptype* buffer, const ievt, const ihel ) <===]
    static constexpr auto ieventAccessIhelConst =
      MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template ieventAccessFieldConst<int>;
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
    // [Signature (non-const, SCALAR) ===> fptype& kernelAccess_s( fptype* buffer, const int ihel ) <===]
    static constexpr auto kernelAccessIhel_s =
      KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessField<int>; // requires cuda 11.4

    // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal)
    // [Signature (non const, SCALAR OR VECTOR) ===> fptype_sv& kernelAccess( const fptype* buffer, const int ihel ) <===]
    static __host__ __device__ inline fptype_sv&
    kernelAccessIhel( fptype* buffer )
    {
      fptype& out = kernelAccessIhel_s( buffer );
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
    // [Signature (const) ===> const fptype& kernelAccessConst( const fptype* buffer, const int ihel ) <===]
    static constexpr auto kernelAccessIhelConst =
      KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessFieldConst<int>; // requires cuda 11.4
  };

  //----------------------------------------------------------------------------

  typedef KernelAccessMatrixElements<false> HostAccessMatrixElements;
  typedef KernelAccessMatrixElements<true> DeviceAccessMatrixElements;

  //----------------------------------------------------------------------------

} // end namespace mg5amcGpu/mg5amcCpu

#endif // MemoryAccessMatrixElements_H
