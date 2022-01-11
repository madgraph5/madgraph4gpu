#ifndef MemoryAccessRandomNumbers_H
#define MemoryAccessRandomNumbers_H 1

#include "mgOnGpuConfig.h"

#include "MemoryAccessHelpers.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for random numbers
// This implementation uses an AOSOA[npagR][nparf][np4][neppR] where nevt=npagR*neppR
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessRandomNumbersBase//_AOSOAv1
{
public:

  // Number of Events Per Page in the random number AOSOA memory buffer layout
  // *** NB Different values of neppR lead to different physics results: the ***
  // *** same 1d array is generated, but it is interpreted in different ways ***
  static constexpr int neppR = 8; // HARDCODED TO GIVE ALWAYS THE SAME PHYSICS RESULTS!
  //static constexpr int neppR = 1; // AOS (tests of sectors/requests)

 private:

  friend class MemoryAccessHelper<MemoryAccessRandomNumbersBase>;
  friend class KernelAccessHelper<MemoryAccessRandomNumbersBase, true>;
  friend class KernelAccessHelper<MemoryAccessRandomNumbersBase, false>;
  
  // The number of components of a 4-momentum
  static constexpr int np4 = mgOnGpu::np4;

  // The number of final state particles in this physics process
  static constexpr int nparf = mgOnGpu::nparf;

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
    constexpr int ip4 = 0;
    constexpr int iparf = 0;
    const int ipagR = ievt/neppR; // #event "R-page"
    const int ieppR = ievt%neppR; // #event in the current event R-page
    return &( buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR] ); // AOSOA[ipagR][iparf][ip4][ieppR]
  }

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // (Non-const memory access to field in an event record)
  static
  __host__ __device__ inline
  fptype& decodeRecord( fptype* buffer,
                        const int ip4,
                        const int iparf )
  {
    constexpr int ipagR = 0;
    constexpr int ieppR = 0;
    return buffer[ipagR*nparf*np4*neppR + iparf*np4*neppR + ip4*neppR + ieppR]; // AOSOA[ipagR][iparf][ip4][ieppR]
  }

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
class MemoryAccessRandomNumbers : public MemoryAccessRandomNumbersBase
{
public:

  // (Non-const memory access to event record from ievent)
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessRandomNumbersBase>::ieventAccessRecord;

  // (Const memory access to event record from ievent)
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessRandomNumbersBase>::ieventAccessRecordConst;

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecordIp4Iparf = MemoryAccessHelper<MemoryAccessRandomNumbersBase>::decodeRecord;

  // [NOTE THE USE OF THE TEMPLATE KEYWORD IN ALL OF THE FOLLOWING TEMPLATE FUNCTION INSTANTIATIONS]
  // (Const memory access to field in an event record)
  static constexpr auto decodeRecordIp4IparfConst =
    MemoryAccessHelper<MemoryAccessRandomNumbersBase>::template decodeRecordConst<int, int>;

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccessIp4Iparf =
    MemoryAccessHelper<MemoryAccessRandomNumbersBase>::template ieventAccessField<int, int>;

  // (Const memory access to field from ievent)
  static constexpr auto ieventAccessIp4IparfConst =
    MemoryAccessHelper<MemoryAccessRandomNumbersBase>::template ieventAccessFieldConst<int, int>;

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
template<bool onDevice>
class KernelAccessRandomNumbers
{
public:

  // (Non-const memory access to field from kernel)
  static constexpr auto kernelAccessIp4Iparf =
    KernelAccessHelper<MemoryAccessRandomNumbersBase, onDevice>::template kernelAccessField<int, int>;

  // (Const memory access to field from kernel)
  static constexpr auto kernelAccessIp4IparfConst =
    KernelAccessHelper<MemoryAccessRandomNumbersBase, onDevice>::template kernelAccessFieldConst<int, int>;

};

//----------------------------------------------------------------------------

typedef KernelAccessRandomNumbers<false> HostAccessRandomNumbers;
typedef KernelAccessRandomNumbers<true> DeviceAccessRandomNumbers;

//----------------------------------------------------------------------------

#endif // MemoryAccessRandomNumbers_H
