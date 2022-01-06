#ifndef MemoryAccessMatrixElements_H
#define MemoryAccessMatrixElements_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
//#include "mgOnGpuVectors.h"

#include "MemoryAccessHelpers.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for matrix elements
// This implementation uses a plain ARRAY[nevt]
// [If many implementations are used, a suffix _ARRAYv1 should be appended to the class name]
class MemoryAccessMatrixElementsBase//_ARRAYv1
{
private:

  friend class MemoryAccessHelper<MemoryAccessMatrixElementsBase>;
  friend class KernelAccessHelper<MemoryAccessMatrixElementsBase, true>;
  friend class KernelAccessHelper<MemoryAccessMatrixElementsBase, false>;

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
    return &( buffer[ievt] ); // ARRAY[nevt]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // (Non-const memory access to field in an event record)
  static
  __host__ __device__ inline
  fptype& decodeRecord( fptype* buffer )
  {
    constexpr int ievt = 0;
    return buffer[ievt]; // ARRAY[nevt]
  }

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
class MemoryAccessMatrixElements : public MemoryAccessMatrixElementsBase
{
public:

  // (Non-const memory access to event record from ievent)
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::ieventAccessRecord;

  // (Const memory access to event record from ievent)
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::ieventAccessRecordConst;

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecord = MemoryAccessHelper<MemoryAccessMatrixElementsBase>::decodeRecord;

  // [NOTE THE USE OF THE TEMPLATE KEYWORD IN ALL OF THE FOLLOWING TEMPLATE FUNCTION INSTANTIATIONS]
  // (Const memory access to field in an event record)
  static constexpr auto decodeRecordConst =
    MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template decodeRecordConst<>;

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccess =
    MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template ieventAccessField<>;

  // (Const memory access to field from ievent)
  static constexpr auto ieventAccessConst =
    MemoryAccessHelper<MemoryAccessMatrixElementsBase>::template ieventAccessFieldConst<>;

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
template<bool onDevice>
class KernelAccessMatrixElements
{
public:

  // (Non-const memory access to field from kernel)
  // FINAL IMPLEMENTATION FOR CUDA 11.4
  //static constexpr auto kernelAccess =
  //  KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessField<>;
  // TEMPORARY HACK FOR CUDA 11.1
  static __host__ __device__ inline
  fptype& kernelAccess( fptype* buffer )
  {
    return KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessField<>( buffer );
  }

  // (Const memory access to field from kernel)
  // FINAL IMPLEMENTATION FOR CUDA 11.4
  //static constexpr auto kernelAccessConst =
  //  KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessFieldConst<>;
  // TEMPORARY HACK FOR CUDA 11.1
  static __host__ __device__ inline
  const fptype& kernelAccessConst( const fptype* buffer )
  {
    return KernelAccessHelper<MemoryAccessMatrixElementsBase, onDevice>::template kernelAccessFieldConst<>( buffer );
  }

};

//----------------------------------------------------------------------------

typedef KernelAccessMatrixElements<false> HostAccessMatrixElements;
typedef KernelAccessMatrixElements<true> DeviceAccessMatrixElements;

//----------------------------------------------------------------------------

#endif // MemoryAccessMatrixElements_H
