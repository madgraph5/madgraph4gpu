//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "mgOnGpuConfig.h"

#include "HelAmps_sm.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

mgDebugDeclare();

#ifdef __CUDACC__
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__
  void FFV1_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    const cxtype_sv TMP0 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) + ( F1[3] * ( F2[4] * ( V3[3]- cI * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) + ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) ) );
    (*vertex) = COUP * - cI * TMP0;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__
  void FFV1P0_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    V3[0] = + F1[0] + F2[0];
    V3[1] = + F1[1] + F2[1];
    const fptype_sv P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype_sv denom = COUP / ( (P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * (- cI) * ( F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3] );
    V3[3] = denom * (- cI) * ( -F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2] );
    V3[4] = denom * (- cI) * ( - cI * ( F1[2] * F2[5] + F1[5] * F2[2] ) + cI * ( F1[3] * F2[4] + F1[4] * F2[3] ) );
    V3[5] = denom * (- cI) * ( -F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2] );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__
  void FFV2_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    const cxtype_sv TMP1 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[3] * ( F2[4] * ( V3[3]- cI * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) );
    (*vertex) = COUP * - cI * TMP1;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__
  void FFV2_3( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = + F1[0] + F2[0];
    V3[1] = + F1[1] + F2[1];
    const fptype_sv P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    const cxtype_sv TMP2 = ( F1[2] * ( F2[4] * ( P3[0] + P3[3] ) + F2[5] * ( P3[1] + cI * ( P3[2] ) ) ) + F1[3] * ( F2[4] * ( P3[1]- cI * ( P3[2] ) ) + F2[5] * ( P3[0] - P3[3] ) ) );
    const cxtype_sv denom = COUP / ( (P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * (- cI) * ( F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2 );
    V3[3] = denom * (- cI) * ( -F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP2 );
    V3[4] = denom * (- cI) * ( - cI * ( F1[2] * F2[5] ) + cI * ( F1[3] * F2[4] ) - P3[2] * OM3 * TMP2 );
    V3[5] = denom * (- cI) * ( -F1[2] * F2[4] - P3[3] * OM3 * TMP2 + F1[3] * F2[5] );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__
  void FFV4_0( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               cxtype_sv* vertex )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    const cxtype_sv TMP1 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[3] * ( F2[4] * ( V3[3]- cI * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) );
    const cxtype_sv TMP3 = ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) );
    (*vertex) = COUP * (- one) * ( +cI * ( TMP1 ) + two * cI * ( TMP3 ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__
  void FFV4_3( const cxtype_sv F1[],
               const cxtype_sv F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype_sv V3[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = + F1[0] + F2[0];
    V3[1] = + F1[1] + F2[1];
    const fptype_sv P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    constexpr fptype two( 2. );
    constexpr fptype half( 1. / 2. );
    const cxtype_sv TMP2 = ( F1[2] * ( F2[4] * ( P3[0] + P3[3] ) + F2[5] * ( P3[1] + cI * ( P3[2] ) ) ) + F1[3] * ( F2[4] * ( P3[1]- cI * ( P3[2] ) ) + F2[5] * ( P3[0] - P3[3] ) ) );
    const cxtype_sv TMP4 = ( F1[4] * ( F2[2] * ( P3[0] - P3[3] ) - F2[3] * ( P3[1] + cI * ( P3[2] ) ) ) + F1[5] * ( F2[2] * ( -P3[1] + cI * ( P3[2] ) ) + F2[3] * ( P3[0] + P3[3] ) ) );
    const cxtype_sv denom = COUP / ( (P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * (- two * cI) * ( OM3 * - half * P3[0] * ( TMP2 + two * ( TMP4 ) ) + ( +half * ( F1[2] * F2[4] + F1[3] * F2[5] ) + F1[4] * F2[2] + F1[5] * F2[3] ) );
    V3[3] = denom * (- two * cI) * ( OM3 * - half * P3[1] * ( TMP2 + two * ( TMP4 ) ) + ( - half * ( F1[2] * F2[5] + F1[3] * F2[4] ) + F1[4] * F2[3] + F1[5] * F2[2] ) );
    V3[4] = denom * two * cI * ( OM3 * half * P3[2] * ( TMP2 + two * ( TMP4 ) ) + ( +half * cI * ( F1[2] * F2[5] )- half * cI * ( F1[3] * F2[4] )- cI * ( F1[4] * F2[3] ) + cI * ( F1[5] * F2[2] ) ) );
    V3[5] = denom * two * cI * ( OM3 * half * P3[3] * ( TMP2 + two * ( TMP4 ) ) + ( +half * ( F1[2] * F2[4] )- half * ( F1[3] * F2[5] ) - F1[4] * F2[2] + F1[5] * F2[3] ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output amplitude 'vertex' from the input wavefunctions F1[6], F2[6], V3[6]
  __device__
  void FFV2_4_0( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype_sv V3[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype_sv* vertex )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    const cxtype_sv TMP1 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[3] * ( F2[4] * ( V3[3]- cI * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) );
    const cxtype_sv TMP3 = ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) );
    (*vertex) = (- one) * ( COUP2 * ( +cI * ( TMP1 ) + two * cI * ( TMP3 ) ) + cI * ( TMP1 * COUP1 ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'V3[6]' from the input wavefunctions F1[6], F2[6]
  __device__
  void FFV2_4_3( const cxtype_sv F1[],
                 const cxtype_sv F2[],
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    const fptype OM3 = ( M3 != 0. ? 1. / ( M3 * M3 ) : 0. );
    V3[0] = + F1[0] + F2[0];
    V3[1] = + F1[1] + F2[1];
    const fptype_sv P3[4] = { -cxreal( V3[0] ), -cxreal( V3[1] ), -cximag( V3[1] ), -cximag( V3[0] ) };
    constexpr fptype one( 1. );
    constexpr fptype two( 2. );
    constexpr fptype half( 1. / 2. );
    const cxtype_sv TMP2 = ( F1[2] * ( F2[4] * ( P3[0] + P3[3] ) + F2[5] * ( P3[1] + cI * ( P3[2] ) ) ) + F1[3] * ( F2[4] * ( P3[1]- cI * ( P3[2] ) ) + F2[5] * ( P3[0] - P3[3] ) ) );
    const cxtype_sv TMP4 = ( F1[4] * ( F2[2] * ( P3[0] - P3[3] ) - F2[3] * ( P3[1] + cI * ( P3[2] ) ) ) + F1[5] * ( F2[2] * ( -P3[1] + cI * ( P3[2] ) ) + F2[3] * ( P3[0] + P3[3] ) ) );
    const cxtype_sv denom = one / ( (P3[0] * P3[0] ) - ( P3[1] * P3[1] ) - ( P3[2] * P3[2] ) - ( P3[3] * P3[3] ) - M3 * ( M3 - cI * W3 ) );
    V3[2] = denom * (- two * cI) * ( COUP2 * ( OM3 * - half * P3[0] * ( TMP2 + two * ( TMP4 ) ) + ( +half * ( F1[2] * F2[4] + F1[3] * F2[5] ) + F1[4] * F2[2] + F1[5] * F2[3] ) ) + half * ( COUP1 * ( F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2 ) ) );
    V3[3] = denom * (- two * cI) * ( COUP2 * ( OM3 * - half * P3[1] * ( TMP2 + two * ( TMP4 ) ) + ( - half * ( F1[2] * F2[5] + F1[3] * F2[4] ) + F1[4] * F2[3] + F1[5] * F2[2] ) )- half * ( COUP1 * ( F1[2] * F2[5] + F1[3] * F2[4] + P3[1] * OM3 * TMP2 ) ) );
    V3[4] = denom * cI * ( COUP2 * ( OM3 * P3[2] * ( TMP2 + two * ( TMP4 ) ) + ( +cI * ( F1[2] * F2[5] )- cI * ( F1[3] * F2[4] )- two * cI * ( F1[4] * F2[3] ) + two * cI * ( F1[5] * F2[2] ) ) ) + COUP1 * ( +cI * ( F1[2] * F2[5] )- cI * ( F1[3] * F2[4] ) + P3[2] * OM3 * TMP2 ) );
    V3[5] = denom * two * cI * ( COUP2 * ( OM3 * half * P3[3] * ( TMP2 + two * ( TMP4 ) ) + ( +half * ( F1[2] * F2[4] )- half * ( F1[3] * F2[5] ) - F1[4] * F2[2] + F1[5] * F2[3] ) ) + half * ( COUP1 * ( F1[2] * F2[4] + P3[3] * OM3 * TMP2 - F1[3] * F2[5] ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

} // end namespace

