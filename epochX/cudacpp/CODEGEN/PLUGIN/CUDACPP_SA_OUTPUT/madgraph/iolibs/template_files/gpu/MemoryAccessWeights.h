#ifndef MemoryAccessWeights_H
#define MemoryAccessWeights_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for weights
// This implementation uses a plain ARRAY[nevt]
// [If many implementations are used, a suffix _ARRAYv1 should be appended to the class name]
class MemoryAccessWeightsBase//_ARRAYv1
{
private:

  friend class MemoryAccessHelper<MemoryAccessWeightsBase>;
  friend class KernelAccessHelper<MemoryAccessWeightsBase, true>;
  friend class KernelAccessHelper<MemoryAccessWeightsBase, false>;

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
class MemoryAccessWeights : public MemoryAccessWeightsBase
{
public:

  // (Non-const memory access to event record from ievent)
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessWeightsBase>::ieventAccessRecord;

  // (Const memory access to event record from ievent)
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessWeightsBase>::ieventAccessRecordConst;

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecord = MemoryAccessHelper<MemoryAccessWeightsBase>::decodeRecord;

  // [NOTE THE USE OF THE TEMPLATE KEYWORD IN ALL OF THE FOLLOWING TEMPLATE FUNCTION INSTANTIATIONS]
  // (Const memory access to field in an event record)
  static constexpr auto decodeRecordConst =
    MemoryAccessHelper<MemoryAccessWeightsBase>::template decodeRecordConst<>;

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccess =
    MemoryAccessHelper<MemoryAccessWeightsBase>::template ieventAccessField<>;

  // (Const memory access to field from ievent)
  static constexpr auto ieventAccessConst =
    MemoryAccessHelper<MemoryAccessWeightsBase>::template ieventAccessFieldConst<>;

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
template<bool onDevice>
class KernelAccessWeights
{
public:

  // (Non-const memory access to field from kernel)
  // FINAL IMPLEMENTATION FOR CUDA 11.4
  //static constexpr auto kernelAccess =
  //  KernelAccessHelper<MemoryAccessWeightsBase, onDevice>::template kernelAccessField<>;
  // TEMPORARY HACK FOR CUDA 11.1
  static __host__ __device__ inline
  fptype& kernelAccess( fptype* buffer )
  {
    return KernelAccessHelper<MemoryAccessWeightsBase, onDevice>::template kernelAccessField<>( buffer );
  }

  // (Const memory access to field from kernel)
  // FINAL IMPLEMENTATION FOR CUDA 11.4
  //static constexpr auto kernelAccessConst =
  //  KernelAccessHelper<MemoryAccessWeightsBase, onDevice>::template kernelAccessFieldConst<>;
  // TEMPORARY HACK FOR CUDA 11.1
  static __host__ __device__ inline
  const fptype& kernelAccessConst( const fptype* buffer )
  {
    return KernelAccessHelper<MemoryAccessWeightsBase, onDevice>::template kernelAccessFieldConst<>( buffer );
  }

};

//----------------------------------------------------------------------------

typedef KernelAccessWeights<false> HostAccessWeights;
typedef KernelAccessWeights<true> DeviceAccessWeights;

//----------------------------------------------------------------------------

#endif // MemoryAccessWeights_H
