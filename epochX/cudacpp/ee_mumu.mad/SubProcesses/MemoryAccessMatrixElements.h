#ifndef MemoryAccessMatrixElements_H
#define MemoryAccessMatrixElements_H 1

#include "mgOnGpuConfig.h"

#include "MemoryAccessHelpers.h"

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

  /*
  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& kernelAccess( fptype* buffer ) <===]
  // FINAL IMPLEMENTATION FOR CUDA 11.4
  static constexpr auto kernelAccess =
    KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessField<>;
  */

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (non-const) ===> fptype& kernelAccess( fptype* buffer ) <===]
  // TEMPORARY HACK FOR CUDA 11.1
  static __host__ __device__ inline fptype&
  kernelAccess( fptype* buffer )
  {
    return KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessField<>( buffer );
  }

  /*
  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const) ===> const fptype& kernelAccessConst( const fptype* buffer ) <===]
  // FINAL IMPLEMENTATION FOR CUDA 11.4
  static constexpr auto kernelAccessConst =
    KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessFieldConst<>;
  */

  // Locate a field (output) in a memory buffer (input) from a kernel event-indexing mechanism (internal) and the given field indexes (input)
  // [Signature (const) ===> const fptype& kernelAccessConst( const fptype* buffer ) <===]
  // TEMPORARY HACK FOR CUDA 11.1
  static __host__ __device__ inline const fptype&
  kernelAccessConst( const fptype* buffer )
  {
    return KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessFieldConst<>( buffer );
  }
};

//----------------------------------------------------------------------------

typedef KernelAccessMatrixElements<false> HostAccessMatrixElements;
typedef KernelAccessMatrixElements<true> DeviceAccessMatrixElements;

//----------------------------------------------------------------------------

#endif // MemoryAccessMatrixElements_H
