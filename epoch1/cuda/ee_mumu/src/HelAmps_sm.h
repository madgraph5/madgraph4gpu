//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H 1

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

namespace MG5_sm
{
  /*
  // Decode momentum AOSOA: compute address of fptype for the given particle, 4-momentum component and event
  // Return the fptype by reference (equivalent to returning its memory address)
  __device__
  inline const fptype& pIparIp4Ievt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event M-page
    //printf( "%2d %2d %8d %8.3f\n", ipar, ip4, ievt, momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
    return momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }
  */

  //--------------------------------------------------------------------------

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
#if defined __CUDACC__ or not defined MGONGPU_CPPSIMD
    __device__ p4type_sv( const p4type_ref ref ) : p0( ref.p0 ), p1( ref.p1 ), p2( ref.p2 ), p3( ref.p3 ){}
#endif
  };

#ifndef __CUDACC__
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
#endif

  // Decode momenta AOSOA: return 4-momentum references of a given particle for a given event
  // Returning by reference seems irrelevant for performance, but allows a simpler code structure
  __device__
  inline p4type_ref p4IparIevt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                const int ipar,
                                const int ievt )
  {
    /*
    return p4type_sv{ pIparIp4Ievt( momenta1d, ipar, 0, ievt ),
                      pIparIp4Ievt( momenta1d, ipar, 1, ievt ),
                      pIparIp4Ievt( momenta1d, ipar, 2, ievt ),
                      pIparIp4Ievt( momenta1d, ipar, 3, ievt ) };
    */
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #event "M-page"
    const int ieppM = ievt%neppM; // #event in the current event "M-page"
    const int index0 = ipagM*npar*np4*neppM + ipar*np4*neppM + ieppM; // 1d-index for AOSOA[ipagM][ipar][0][ieppM]
    //for ( ip4=0; ip4<np4; ip4++) printf( "%2d %2d %8d %8.3f\n", ipar, ip4, ievt, momenta1d[index0 + ip4*neppM] );
    return p4type_ref{ momenta1d[index0 + 0*neppM],   // AOSOA[ipagM][ipar][0][ieppM]
                       momenta1d[index0 + 1*neppM],   // AOSOA[ipagM][ipar][1][ieppM]
                       momenta1d[index0 + 2*neppM],   // AOSOA[ipagM][ipar][2][ieppM]
                       momenta1d[index0 + 3*neppM] }; // AOSOA[ipagM][ipar][3][ieppM]
  }

#ifndef __CUDACC__
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
    // Use c++17 "if constexpr": compile-time branching
    if constexpr ( ( neppM >= neppV ) && ( neppM%neppV == 0 ) )
    {
      constexpr bool useReinterpretCastIfPossible = true; // DEFAULT
      //constexpr bool useReinterpretCastIfPossible = false; // FOR PERFORMANCE TESTS
      constexpr bool skipAlignmentCheck = true; // DEFAULT (MAY SEGFAULT, NEEDS A SANITY CHECK ELSEWHERE!)
      //constexpr bool skipAlignmentCheck = false; // SLOWER BUT SAFER
      if constexpr ( useReinterpretCastIfPossible && skipAlignmentCheck )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! skip alignment check" << std::endl; first=false; } // SLOWS DOWN...
        // Fastest (4.92E6 in eemumu 512y)
        // This assumes alignment for momenta1d without checking - causes segmentation fault in reinterpret_cast if not aligned!
        return p4typevFromAlignedArray( p4IparIevt( momenta1d, ipar, ievt0 ) ); // use reinterpret_cast
      }
      else if ( useReinterpretCastIfPossible && ( (size_t)(momenta1d) % mgOnGpu::cppAlign == 0 ) )
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! alignment ok, use reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A bit (6%) slower (4.62E6 in eemumu 512y) because of the alignment check
        // This explicitly checks alignment for momenta1d to avoid segmentation faults in reinterpret_cast
        return p4typevFromAlignedArray( p4IparIevt( momenta1d, ipar, ievt0 ) ); // use reinterpret_cast
      }
      else
      {
        //static bool first=true; if( first ){ std::cout << "WARNING! AOSOA but no reinterpret cast" << std::endl; first=false; } // SLOWS DOWN...
        // A bit (2%) slower (4.86E6 in eemumu 512y)
        // This does not require alignment for momenta1d, but it requires AOSOA with neppM>=neppV and neppM%neppV==0
        return p4typevFromUnalignedArray( p4IparIevt( momenta1d, ipar, ievt0 ) ); // do not use reinterpret_cast
      }
    }
    else
    {
      // Much (20%) slower (4.07E6 in eemumu 512y)
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

#ifdef MGONGPU_INLINE_HELAMPS
#define INLINE inline
#define ALWAYS_INLINE __attribute__((always_inline))
#else
#define INLINE
#define ALWAYS_INLINE
#endif

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixxxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ipzxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void imzxxx( const p4type_sv& p4vec,        // input: four-momentum
               //const fptype fmass,          // ASSUME fmass==0
               const int nhel,                // input: -1 or +1 (helicity of fermion)
               const int nsf,                 // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi ) ALWAYS_INLINE; // output: wavefunction[(nw6==6)]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void ixzxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fi,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void vxxxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               const fptype vmass,
               const int nhel,                 // input: -1, 0 (only if vmass!=0) or +1 (helicity of vector boson)
               const int nsv,                  // input: +1 (final) or -1 (initial)
               cxtype_sv* vc,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void sxxxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               const fptype,                   // WARNING: "smass" unused (missing in Fortran)
               const int,                      // WARNING: "nhel" unused (missing in Fortran) - scalar has no helicity
               const int nss,                  // input: +1 (final) or -1 (initial)
               cxtype_sv sc[3],                // output: wavefunction[3] - not [6], this is for scalars
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxxxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               const fptype fmass,
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void opzxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void omzxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void oxzxxx( const fptype* allmomenta,       // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,           // ASSUME fmass==0
               const int nhel,                 // input: -1 or +1 (helicity of fermion)
               const int nsf,                  // input: +1 (particle) or -1 (antiparticle)
               cxtype_sv* fo,                  // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar ) ALWAYS_INLINE; // input: particle# out of npar

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV1_0( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv V3[],                // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;   // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV1P0_3( const cxtype_sv F1[],              // input: wavefunction1[6]
                 const cxtype_sv F2[],              // input: wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;    // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_0( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv V3[],                // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;   // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_3( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] ) ALWAYS_INLINE;      // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV4_0( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv V3[],                // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex ) ALWAYS_INLINE;   // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV4_3( const cxtype_sv F1[],                // input: wavefunction1[6]
               const cxtype_sv F2[],                // input: wavefunction2[6]
               const cxtype_sv COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] ) ALWAYS_INLINE;      // output: wavefunction3[6]

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_4_0( const cxtype_sv F1[],              // input: wavefunction1[6]
                 const cxtype_sv F2[],              // input: wavefunction2[6]
                 const cxtype_sv V3[],              // input: wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype_sv* vertex ) ALWAYS_INLINE; // output: amplitude

  //--------------------------------------------------------------------------

  __device__ INLINE
  void FFV2_4_3( const cxtype_sv F1[],              // input: wavefunction1[6]
                 const cxtype_sv F2[],              // input: wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] ) ALWAYS_INLINE;    // output: wavefunction3[6]

  //--------------------------------------------------------------------------

} // end namespace

#endif // HelAmps_sm_H
