#ifndef MemoryAccessMomenta_H
#define MemoryAccessMomenta_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include "MemoryAccess.h" // FIXME! rename as MemoryAccessVectors.h eventually
#include "MemoryAccessHelpers.h"

//----------------------------------------------------------------------------

// A class describing the internal layout of memory buffers for momenta
// This implementation uses an AOSOA[npagM][npar][np4][neppM] where nevt=npagM*neppM
// [If many implementations are used, a suffix _AOSOAv1 should be appended to the class name]
class MemoryAccessMomentaBase//_AOSOAv1
{
public:

  // Number of Events Per Page in the momenta AOSOA memory buffer layout
  static constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time

private:

  friend class MemoryAccessHelper<MemoryAccessMomentaBase>;
  friend class KernelAccessHelper<MemoryAccessMomentaBase, true>;
  friend class KernelAccessHelper<MemoryAccessMomentaBase, false>;
  
  // The number of components of a 4-momentum
  static constexpr int np4 = mgOnGpu::np4;

  // The number of particles in this physics process
  static constexpr int npar = mgOnGpu::npar;

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
    constexpr int ipar = 0;
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    return &( buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] ); // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  //--------------------------------------------------------------------------

  // Locate a field (output) of an event record (input) from the given field indexes (input)
  // (Non-const memory access to field in an event record)
  static
  __host__ __device__ inline
  fptype& decodeRecord( fptype* buffer,
                        const int ip4,
                        const int ipar )
  {
    constexpr int ipagM = 0;
    constexpr int ieppM = 0;
    return buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on explicit event numbers
class MemoryAccessMomenta : public MemoryAccessMomentaBase
{
public:

  // (Non-const memory access to event record from ievent)
  static constexpr auto ieventAccessRecord = MemoryAccessHelper<MemoryAccessMomentaBase>::ieventAccessRecord;

  // (Const memory access to event record from ievent)
  static constexpr auto ieventAccessRecordConst = MemoryAccessHelper<MemoryAccessMomentaBase>::ieventAccessRecordConst;

  // (Non-const memory access to field in an event record)
  static constexpr auto decodeRecordIp4Ipar = MemoryAccessHelper<MemoryAccessMomentaBase>::decodeRecord;

  // [NOTE THE USE OF THE TEMPLATE KEYWORD IN ALL OF THE FOLLOWING TEMPLATE FUNCTION INSTANTIATIONS]
  // (Const memory access to field in an event record)
  static constexpr auto decodeRecordIp4IparConst =
    MemoryAccessHelper<MemoryAccessMomentaBase>::template decodeRecordConst<int, int>;

  // (Non-const memory access to field from ievent)
  static constexpr auto ieventAccessIp4Ipar =
    MemoryAccessHelper<MemoryAccessMomentaBase>::template ieventAccessField<int, int>;

  // (Const memory access to field from ievent)
  static constexpr auto ieventAccessIp4IparConst =
    MemoryAccessHelper<MemoryAccessMomentaBase>::template ieventAccessFieldConst<int, int>;
  /*
  // (Const memory access to field from ievent - DEBUG version with printouts)
  static
  __host__ __device__ inline
  const fptype& ieventAccessIp4IparConst( const fptype* buffer,
                                          const int ievt,
                                          const int ip4,
                                          const int ipar )
  {
    const fptype& out = MemoryAccessHelper<MemoryAccessMomentaBase>::template ieventAccessFieldConst<int, int>( buffer, ievt, ip4, ipar );
    printf( "ipar=%2d ip4=%2d ievt=%8d out=%8.3f\n", ipar, ip4, ievt, out );
    return out;
  }
  */

};

//----------------------------------------------------------------------------

// A class providing access to memory buffers for a given event, based on implicit kernel rules
template<bool onDevice>
class KernelAccessMomenta
{
public:

  // (Non-const memory access to field from kernel)
  static constexpr auto kernelAccessIp4Ipar =
    KernelAccessHelper<MemoryAccessMomentaBase, onDevice>::template kernelAccessField<int, int>;

  // (Const memory access to field from kernel, scalar)
  static constexpr auto kernelAccessIp4IparConst_s =
    KernelAccessHelper<MemoryAccessMomentaBase, onDevice>::template kernelAccessFieldConst<int, int>;
  /*
  // (Const memory access to field from kernel, scalar - DEBUG version with printouts)
  static
  __host__ __device__ inline
  const fptype& kernelAccessIp4IparConst_s( const fptype* buffer,
                                            const int ip4,
                                            const int ipar )
  {
    const fptype& out = KernelAccessHelper<MemoryAccessMomentaBase, onDevice>::template kernelAccessFieldConst<int, int>( buffer, ip4, ipar );
    printf( "ipar=%2d ip4=%2d ievt=  kernel out=%8.3f\n", ipar, ip4, out );
    return out;
  }
  */

  // (Const memory access to field from kernel, scalar or vector)
  // [FIXME? Eventually return by reference and support aligned arrays only?]
  // [Currently return by value to support also unaligned and arbitrary arrays]
  static
  __host__ __device__ inline
  fptype_sv kernelAccessIp4IparConst( const fptype* buffer,
                                      const int ip4,
                                      const int ipar )
  {
    const fptype& out = kernelAccessIp4IparConst_s( buffer, ip4, ipar );
#ifndef MGONGPU_CPPSIMD
    return out;
#else
    constexpr int neppM = MemoryAccessMomentaBase::neppM;
    using namespace MG5_sm;
    constexpr bool useContiguousEventsIfPossible = true; // DEFAULT
    //constexpr bool useContiguousEventsIfPossible = false; // FOR PERFORMANCE TESTS (treat as arbitrary array even if it is an AOSOA)
    // Use c++17 "if constexpr": compile-time branching
    if constexpr ( useContiguousEventsIfPossible && ( neppM >= neppV ) && ( neppM%neppV == 0 ) )
    {
      //constexpr bool skipAlignmentCheck = true; // FASTEST (MAY SEGFAULT, NEEDS A SANITY CHECK ELSEWHERE!)
      constexpr bool skipAlignmentCheck = false; // NEW DEFAULT: A BIT SLOWER BUT SAFER [UNCOMMENT OUT TO TEST MISALIGNED ACCESS]
      //if constexpr ( skipAlignmentCheck )
      if ( skipAlignmentCheck )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! assume aligned AOSOA, skip check" << std::endl; first=false; } // SLOWS DOWN...
        // Fastest (4.85E6 in eemumu 512y - was 4.93E6 without kernelAccess functions)
        // This assumes alignment for momenta1d without checking - causes segmentation fault in reinterpret_cast if not aligned!
        return fptypevFromAlignedArray( out ); // use reinterpret_cast
      }
      else if ( (size_t)(buffer) % mgOnGpu::cppAlign == 0 )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! aligned AOSOA, reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A tiny bit (<1%) slower because of the alignment check (4.83E6 in eemumu 512y)
        // This explicitly checks buffer alignment to avoid segmentation faults in reinterpret_cast
        return fptypevFromAlignedArray( out ); // use reinterpret_cast
      }
      else
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! AOSOA but no reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A bit (3%) slower (4.70E6 in eemumu 512y)
        // This does not require buffer alignment, but it requires AOSOA with neppM>=neppV and neppM%neppV==0
        return fptypevFromUnalignedArray( out ); // do not use reinterpret_cast
      }
    }
    else
    {
      //static bool first=true; if( first ){ std::cout << "WARNING! arbitrary array" << std::endl; first=false; } // SLOWS DOWN...
      // Much (7-12%) slower (4.30E6 for AOSOA, 4.53E6 for AOS in eemumu 512y)
      //... to do? implementation based on fptypevFromArbitraryArray ...
      std::cout << "ERROR! useContiguousEventsIfPossible=" << useContiguousEventsIfPossible
                << ", neppM=" << neppM << ", neppV=" << neppV << std::endl;
      throw std::logic_error( "MemoryAccessMomenta requires an AOSOA and does not support arbitrary arrays" ); // no path to this statement
    }
#endif
  }

};

//----------------------------------------------------------------------------

typedef KernelAccessMomenta<false> HostAccessMomenta;
typedef KernelAccessMomenta<true> DeviceAccessMomenta;

//----------------------------------------------------------------------------

#endif // MemoryAccessMomenta_H
