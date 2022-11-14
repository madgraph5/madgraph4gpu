//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef MemoryAccess_H
#define MemoryAccess_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

namespace MG5_sm
{
  // =============================================================================
  // *** Generic pattern: indexingFunction( buffer, ievt, additional_indexes ) ***
  // =============================================================================

  // Memory addressing/indexing function (WITH an explicit event number) for the array/buffer that contains momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the 4-momenta for one specific event, given its event number
  // (Non-const memory access)
  __device__ inline
  fptype& indexingFunctionMomenta( fptype* buffer,
                                   const int ievt,
                                   const int ip4
                                   , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
                                   )
  {
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    //printf( "%2d %2d %8d %8.3f\n", ipar, 0, ievt, buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
    return buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  // (Const memory access)
  __device__ inline
  const fptype& indexingFunctionConstMomenta( const fptype* buffer,
                                              const int ievt,
                                              const int ip4
                                              , const int ipar
                                             )
  {
    return indexingFunctionMomenta( const_cast<fptype*>( buffer ), ievt, ip4, ipar );
  }

  //--------------------------------------------------------------------------

#ifndef __CUDACC__
  // Four-momentum references of one particle, for one event
  // All references point to the original buffer of momenta for all particles in all events
  struct p4type_ref
  {
    const fptype& p0;
    const fptype& p1;
    const fptype& p2;
    const fptype& p3;
  };

  // Four-momentum values of one particle, for one event or for one "V-page" of neppV events
  struct p4type_sv
  {
    fptype_sv p0;
    fptype_sv p1;
    fptype_sv p2;
    fptype_sv p3;
#ifndef MGONGPU_CPPSIMD
    __device__ p4type_sv( const p4type_ref ref ) : p0( ref.p0 ), p1( ref.p1 ), p2( ref.p2 ), p3( ref.p3 ){}
#endif
    operator const fptype_sv*() const{ return &p0; }
    operator fptype_sv*(){ return &p0; }
  };

#ifdef MGONGPU_CPPSIMD
  // Build four fptype_v (four vectors of neppV fptype values) from four fptype references,
  // assuming that "pointer(evt#0)+1" indicates "pointer(evt#1)", and that the arrays are aligned
  inline p4type_sv p4typevFromAlignedArray( const p4type_ref ref )
  {
    return p4type_sv{ *reinterpret_cast<const fptype_sv*>( &( ref.p0 ) ),
                      *reinterpret_cast<const fptype_sv*>( &( ref.p1 ) ),
                      *reinterpret_cast<const fptype_sv*>( &( ref.p2 ) ),
                      *reinterpret_cast<const fptype_sv*>( &( ref.p3 ) ) };
  }

  // Build one fptype_v (one vector of neppV fptype values) from one fptype references,
  // assuming that "pointer(evt#0)+1" indicates "pointer(evt#1)", but that the arrays are not aligned
  inline fptype_v fptypevFromUnalignedArray( const fptype& ref )
  {
#if MGONGPU_CPPSIMD == 2
    return fptype_v{ *( &ref   ),
                     *( &ref+1 ) };
#elif MGONGPU_CPPSIMD == 4
    return fptype_v{ *( &ref   ),
                     *( &ref+1 ),
                     *( &ref+2 ),
                     *( &ref+3 ) };
#elif MGONGPU_CPPSIMD == 8
    return fptype_v{ *( &ref   ),
                     *( &ref+1 ),
                     *( &ref+2 ),
                     *( &ref+3 ),
                     *( &ref+4 ),
                     *( &ref+5 ),
                     *( &ref+6 ),
                     *( &ref+7 ) };
#elif MGONGPU_CPPSIMD == 16
    return fptype_v{ *( &ref   ),
                     *( &ref+1 ),
                     *( &ref+2 ),
                     *( &ref+3 ),
                     *( &ref+4 ),
                     *( &ref+5 ),
                     *( &ref+6 ),
                     *( &ref+7 ),
                     *( &ref+8 ),
                     *( &ref+9 ),
                     *( &ref+10 ),
                     *( &ref+11 ),
                     *( &ref+12 ),
                     *( &ref+13 ),
                     *( &ref+14 ),
                     *( &ref+15 ) };
#else
#error Internal error! Unknown MGONGPU_CPPSIMD value
#endif
  }
  
  // Build four fptype_v (four vectors of neppV fptype values) from four fptype references,
  // assuming that "reference(evt#0)+1" indicates "reference(evt#1)", but that the arrays are not aligned
  inline p4type_sv p4typevFromUnalignedArray( const p4type_ref ref )
  {
    return p4type_sv{ fptypevFromUnalignedArray( ref.p0 ),
                      fptypevFromUnalignedArray( ref.p1 ),
                      fptypevFromUnalignedArray( ref.p2 ),
                      fptypevFromUnalignedArray( ref.p3 ) };
  }

  // Build one fptype_v (one vector of neppV fptype values) from one fptype references,
  // with no a priori assumption on how the input fptype array should be decoded
  template <typename Functor>
  inline fptype_v fptypevFromArbitraryArray( Functor decoderIeppv )
  {
#if MGONGPU_CPPSIMD == 2
    return fptype_v{ decoderIeppv( 0 ),
                     decoderIeppv( 1 ) };
#elif MGONGPU_CPPSIMD == 4
    return fptype_v{ decoderIeppv( 0 ),
                     decoderIeppv( 1 ),
                     decoderIeppv( 2 ),
                     decoderIeppv( 3 ) };
#elif MGONGPU_CPPSIMD == 8
    return fptype_v{ decoderIeppv( 0 ),
                     decoderIeppv( 1 ),
                     decoderIeppv( 2 ),
                     decoderIeppv( 3 ),
                     decoderIeppv( 4 ),
                     decoderIeppv( 5 ),
                     decoderIeppv( 6 ),
                     decoderIeppv( 7 ) };
#elif MGONGPU_CPPSIMD == 16
    return fptype_v{ decoderIeppv( 0 ),
                     decoderIeppv( 1 ),
                     decoderIeppv( 2 ),
                     decoderIeppv( 3 ),
                     decoderIeppv( 4 ),
                     decoderIeppv( 5 ),
                     decoderIeppv( 6 ),
                     decoderIeppv( 7 ),
                     decoderIeppv( 8 ),
                     decoderIeppv( 9 ),
                     decoderIeppv( 10 ),
                     decoderIeppv( 11 ),
                     decoderIeppv( 12 ),
                     decoderIeppv( 13 ),
                     decoderIeppv( 14 ),
                     decoderIeppv( 15 ) };
#else
#error Internal error! Unknown MGONGPU_CPPSIMD value
#endif
  }

  // Build four fptype_v (four vectors of neppV fptype values) from four fptype references,
  // with no a priori assumption on how the input fptype array should be decoded
  template <typename Functor0, typename Functor1, typename Functor2, typename Functor3>
  inline p4type_sv p4typevFromArbitraryArray( Functor0 p0decoderIeppv,
                                              Functor1 p1decoderIeppv,
                                              Functor2 p2decoderIeppv,
                                              Functor3 p3decoderIeppv )
  {
    return p4type_sv{ fptypevFromArbitraryArray( p0decoderIeppv ),
                      fptypevFromArbitraryArray( p1decoderIeppv ),
                      fptypevFromArbitraryArray( p2decoderIeppv ),
                      fptypevFromArbitraryArray( p3decoderIeppv ) };
  }
#endif

  // Decode momenta AOSOA: return 4-momentum references of a given particle for a given event
  // Returning by reference seems irrelevant for performance, but allows a simpler code structure
  __device__
  inline p4type_ref p4IparIevt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                const int ipar,
                                const int ievt )
  {
    return p4type_ref{ indexingFunctionConstMomenta( momenta1d, ievt, 0, ipar ),
                       indexingFunctionConstMomenta( momenta1d, ievt, 1, ipar ),
                       indexingFunctionConstMomenta( momenta1d, ievt, 2, ipar ),
                       indexingFunctionConstMomenta( momenta1d, ievt, 3, ipar ) };
  }

  // Decode momenta AOSOA: return 4-momentum of a given particle for a given event (or 4-momenta for one "V-page" of neppV events)
  // Return by value one fptype (or one fptype_v SIMD vector of neppV fptype values) for each 4-momentum component
  // Strictly speaking, returning the fptype_v by value is only unavoidable for neppM<neppV
  // For neppM>=neppV (both being powers of 2), the momentum neppM-AOSOA is reinterpreted in terms of neppV-vectors:
  // it could also be returned by reference, but no performance degradation is observed when returning by value
  inline p4type_sv p4IparIpagV( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                const int ipar,
                                const int ipagV )
  {
    /*
    return p4type_sv{ pIparIp4IpagV( momenta1d, ipar, 0, ipagV ),
                      pIparIp4IpagV( momenta1d, ipar, 1, ipagV ),
                      pIparIp4IpagV( momenta1d, ipar, 2, ipagV ),
                      pIparIp4IpagV( momenta1d, ipar, 3, ipagV ) };
    */
    const int ievt0 = ipagV*neppV; // virtual event V-page ipagV contains neppV events [ievt0...ievt0+neppV-1]
#ifdef MGONGPU_CPPSIMD
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    constexpr bool useContiguousEventsIfPossible = true; // DEFAULT
    //constexpr bool useContiguousEventsIfPossible = false; // FOR PERFORMANCE TESTS (treat as arbitrary array even if it is an AOSOA)
    // Use c++17 "if constexpr": compile-time branching
    if constexpr ( useContiguousEventsIfPossible && ( neppM >= neppV ) && ( neppM%neppV == 0 ) )
    {
      //constexpr bool skipAlignmentCheck = true; // FASTEST (MAY SEGFAULT, NEEDS A SANITY CHECK ELSEWHERE!)
      constexpr bool skipAlignmentCheck = false; // NEW DEFAULT: A BIT SLOWER BUT SAFER [UNCOMMENT OUT TO TEST MISALIGNED ACCESS]
      if constexpr ( skipAlignmentCheck )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! skip alignment check" << std::endl; first=false; } // SLOWS DOWN...
        // Fastest (4.85E6 in eemumu 512y - was 4.93E6 without kernelAccess functions)
        // This assumes alignment for momenta1d without checking - causes segmentation fault in reinterpret_cast if not aligned!
        return p4typevFromAlignedArray( p4IparIevt( momenta1d, ipar, ievt0 ) ); // use reinterpret_cast
      }
      else if ( (size_t)(momenta1d) % mgOnGpu::cppAlign == 0 )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! alignment ok, use reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A tiny bit (<1%) slower because of the alignment check (4.83E6 in eemumu 512y)
        // This explicitly checks alignment for momenta1d to avoid segmentation faults in reinterpret_cast
        return p4typevFromAlignedArray( p4IparIevt( momenta1d, ipar, ievt0 ) ); // use reinterpret_cast
      }
      else
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! AOSOA but no reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A bit (3%) slower (4.70E6 in eemumu 512y)
        // This does not require alignment for momenta1d, but it requires AOSOA with neppM>=neppV and neppM%neppV==0
        return p4typevFromUnalignedArray( p4IparIevt( momenta1d, ipar, ievt0 ) ); // do not use reinterpret_cast
      }
    }
    else
    {
      static bool first=true; if( first ){ std::cout << "WARNING! arbitrary array" << std::endl; first=false; } // SLOWS DOWN...
      // Much (7-12%) slower (4.30E6 for AOSOA, 4.53E6 for AOS in eemumu 512y)
      // This does not even require AOSOA with neppM>=neppV and neppM%neppV==0 (e.g. can be used with AOS neppM==1)
      auto p0decoderIeppv = [momenta1d, ipar, ievt0](int ieppV) -> const fptype& { return p4IparIevt( momenta1d, ipar, ievt0+ieppV ).p0; };
      auto p1decoderIeppv = [momenta1d, ipar, ievt0](int ieppV) -> const fptype& { return p4IparIevt( momenta1d, ipar, ievt0+ieppV ).p1; };
      auto p2decoderIeppv = [momenta1d, ipar, ievt0](int ieppV) -> const fptype& { return p4IparIevt( momenta1d, ipar, ievt0+ieppV ).p2; };
      auto p3decoderIeppv = [momenta1d, ipar, ievt0](int ieppV) -> const fptype& { return p4IparIevt( momenta1d, ipar, ievt0+ieppV ).p3; };
      return p4typevFromArbitraryArray( p0decoderIeppv, p1decoderIeppv, p2decoderIeppv, p3decoderIeppv );
    }
#else
    return p4IparIevt( momenta1d, ipar, ievt0 );
#endif
  }
#endif

  //--------------------------------------------------------------------------

  // =============================================================================
  // *** Generic pattern: kernelAccessFunction( buffer, additional_indexes ) ***
  // =============================================================================

  // Kernel access function (WITHOUT an explicit event number) for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the 4-momenta for one event or one SIMD vector of events
  // (Non-const memory access)
  __device__ inline
  fptype_sv& kernelAccessMomenta( fptype_sv* buffer,
                                  const int ip4
#ifdef __CUDACC__
                                  , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
#endif
                                  )
  {
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
    //printf( "kernelAccessMomenta: ievt=%d threadId=%d\n", ievt, threadIdx.x );
    return indexingFunctionMomenta( buffer, ievt, ip4, ipar ); // NB fptype and fptype_sv coincide for CUDA
#else
    return buffer[ip4];
#endif
  }

  // (Const memory access)
  __device__ inline
  const fptype_sv& kernelAccessConstMomenta( const fptype_sv* buffer,
                                             const int ip4
#ifdef __CUDACC__
                                             , const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
#endif
                                             )
  {
#ifdef __CUDACC__
    return kernelAccessMomenta( const_cast<fptype_sv*>( buffer ), ip4, ipar );
#else
    return kernelAccessMomenta( const_cast<fptype_sv*>( buffer ), ip4 );
#endif
  }

  //--------------------------------------------------------------------------

  // Kernel access function (WITHOUT an explicit event number) for amplitudes
  // Input: a memory buffer for an arbitrary number of events
  // Output: the amplitude for one event or one SIMD vector of events
  // (Non-const memory access)
  __device__ inline
  cxtype_sv& kernelAccessAmplitudes( cxtype_sv* buffer )
  {
    return *buffer;
  }

  //--------------------------------------------------------------------------

}

#endif // MemoryAccess_H
