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
  SYCL_EXTERNAL inline
  fptype& indexingFunctionMomenta(
          fptype* buffer,
          const int ievt,
          const int ip4,
          const int ipar // TEMPORARY? Move to SOAOSOA? (#309)
          ) {
      using mgOnGpu::np4;
      using mgOnGpu::npar;
      constexpr int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
      const int ipagM = ievt/neppM; // #event "M-page"
      const int ieppM = ievt%neppM; // #event in the current event M-page
      //printf( "%2d %2d %8d %8.3f\n", ipar, 0, ievt, buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
      return buffer[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }

  // (Const memory access)
  SYCL_EXTERNAL inline
  const fptype& indexingFunctionConstMomenta(
          const fptype* buffer,
          const int ievt,
          const int ip4,
          const int ipar
          ) {
      return indexingFunctionMomenta( const_cast<fptype*>( buffer ), ievt, ip4, ipar );
  }

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

    SYCL_EXTERNAL p4type_sv( const p4type_ref ref ) : p0( ref.p0 ), p1( ref.p1 ), p2( ref.p2 ), p3( ref.p3 ){}
    operator const fptype_sv*() const{ return &p0; }
    operator fptype_sv*(){ return &p0; }
  };

  // Decode momenta AOSOA: return 4-momentum references of a given particle for a given event
  // Returning by reference seems irrelevant for performance, but allows a simpler code structure
  SYCL_EXTERNAL
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
    return p4IparIevt( momenta1d, ipar, ievt0 );
  }

  //--------------------------------------------------------------------------

  // =============================================================================
  // *** Generic pattern: kernelAccessFunction( buffer, additional_indexes ) ***
  // =============================================================================

  // Kernel access function (WITHOUT an explicit event number) for momenta
  // Input: a memory buffer for an arbitrary number of events
  // Output: the 4-momenta for one event or one SIMD vector of events
  // (Non-const memory access)
  SYCL_EXTERNAL inline
  fptype_sv& kernelAccessMomenta( fptype_sv* buffer,
                                  const int ip4
                                  ) {
    return buffer[ip4];
  }

  // (Const memory access)
  SYCL_EXTERNAL inline
  const fptype_sv& kernelAccessConstMomenta( const fptype_sv* buffer,
                                             const int ip4
                                             ) {
    return kernelAccessMomenta( const_cast<fptype_sv*>( buffer ), ip4 );
  }

  //--------------------------------------------------------------------------

}

#endif // MemoryAccess_H
