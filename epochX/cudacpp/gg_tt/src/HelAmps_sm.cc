//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 2.9.5, 2021-08-22
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"
#include "HelAmps_sm.h"

//#include <cassert>
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

  // Compute the output wavefunction 'V1[6]' from the input wavefunctions V2[6], V3[6]
  __device__
  void VVV1P0_1( const cxtype_sv V2[],
                 const cxtype_sv V3[],
                 const cxtype COUP,
                 const fptype M1,
                 const fptype W1,
                 cxtype_sv V1[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    const fptype_sv P2[4] = { +cxreal( V2[0] ), +cxreal( V2[1] ), +cximag( V2[1] ), +cximag( V2[0] ) };
    const fptype_sv P3[4] = { +cxreal( V3[0] ), +cxreal( V3[1] ), +cximag( V3[1] ), +cximag( V3[0] ) };
    V1[0] = + V2[0] + V3[0];
    V1[1] = + V2[1] + V3[1];
    const fptype_sv P1[4] = { -cxreal( V1[0] ), -cxreal( V1[1] ), -cximag( V1[1] ), -cximag( V1[0] ) };
    const cxtype_sv TMP0 = ( V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3] );
    const cxtype_sv TMP1 = ( V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3] );
    const cxtype_sv TMP2 = ( P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5] );
    const cxtype_sv TMP3 = ( V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3] );
    const cxtype_sv TMP4 = ( V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5] );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    V1[2] = denom * ( TMP4 * ( - cI * ( P2[0] ) + cI * ( P3[0] ) ) + ( V2[2] * ( - cI * ( TMP0 ) + cI * ( TMP1 ) ) + V3[2] * ( +cI * ( TMP2 )- cI * ( TMP3 ) ) ) );
    V1[3] = denom * ( TMP4 * ( - cI * ( P2[1] ) + cI * ( P3[1] ) ) + ( V2[3] * ( - cI * ( TMP0 ) + cI * ( TMP1 ) ) + V3[3] * ( +cI * ( TMP2 )- cI * ( TMP3 ) ) ) );
    V1[4] = denom * ( TMP4 * ( - cI * ( P2[2] ) + cI * ( P3[2] ) ) + ( V2[4] * ( - cI * ( TMP0 ) + cI * ( TMP1 ) ) + V3[4] * ( +cI * ( TMP2 )- cI * ( TMP3 ) ) ) );
    V1[5] = denom * ( TMP4 * ( - cI * ( P2[3] ) + cI * ( P3[3] ) ) + ( V2[5] * ( - cI * ( TMP0 ) + cI * ( TMP1 ) ) + V3[5] * ( +cI * ( TMP2 )- cI * ( TMP3 ) ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

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
    const cxtype_sv TMP5 = ( F1[2] * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) + ( F1[3] * ( F2[4] * ( V3[3]- cI * ( V3[4] ) ) + F2[5] * ( V3[2] - V3[5] ) ) + ( F1[4] * ( F2[2] * ( V3[2] - V3[5] ) - F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) + F1[5] * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) ) );
    (*vertex) = COUP * - cI * TMP5;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F1[6]' from the input wavefunctions F2[6], V3[6]
  __device__
  void FFV1_1( const cxtype_sv F2[],
               const cxtype_sv V3[],
               const cxtype COUP,
               const fptype M1,
               const fptype W1,
               cxtype_sv F1[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    F1[0] = + F2[0] + V3[0];
    F1[1] = + F2[1] + V3[1];
    const fptype_sv P1[4] = { -cxreal( F1[0] ), -cxreal( F1[1] ), -cximag( F1[1] ), -cximag( F1[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( (P1[0] * P1[0] ) - ( P1[1] * P1[1] ) - ( P1[2] * P1[2] ) - ( P1[3] * P1[3] ) - M1 * ( M1 - cI * W1 ) );
    F1[2] = denom * cI * ( F2[2] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3]- cI * ( V3[4] ) ) + ( P1[2] * ( +cI * ( V3[3] ) + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[3] + cI * ( V3[4] ) ) + ( P1[1] * (- one) * ( V3[2] + V3[5] ) + ( P1[2] * (- one) * ( +cI * ( V3[2] + V3[5] ) ) + P1[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + M1 * ( F2[4] * ( V3[2] + V3[5] ) + F2[5] * ( V3[3] + cI * ( V3[4] ) ) ) ) );
    F1[3] = denom * (- cI) * ( F2[2] * ( P1[0] * ( -V3[3] + cI * ( V3[4] ) ) + ( P1[1] * ( V3[2] - V3[5] ) + ( P1[2] * ( - cI * ( V3[2] ) + cI * ( V3[5] ) ) + P1[3] * ( V3[3]- cI * ( V3[4] ) ) ) ) ) + ( F2[3] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P1[2] * ( +cI * ( V3[3] ) - V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + M1 * ( F2[4] * ( -V3[3] + cI * ( V3[4] ) ) + F2[5] * ( -V3[2] + V3[5] ) ) ) );
    F1[4] = denom * (- cI) * ( F2[4] * ( P1[0] * ( V3[2] + V3[5] ) + ( P1[1] * ( -V3[3] + cI * ( V3[4] ) ) + ( P1[2] * (- one) * ( +cI * ( V3[3] ) + V3[4] ) - P1[3] * ( V3[2] + V3[5] ) ) ) ) + ( F2[5] * ( P1[0] * ( V3[3] + cI * ( V3[4] ) ) + ( P1[1] * ( -V3[2] + V3[5] ) + ( P1[2] * ( - cI * ( V3[2] ) + cI * ( V3[5] ) ) - P1[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + M1 * ( F2[2] * ( -V3[2] + V3[5] ) + F2[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) );
    F1[5] = denom * cI * ( F2[4] * ( P1[0] * ( -V3[3] + cI * ( V3[4] ) ) + ( P1[1] * ( V3[2] + V3[5] ) + ( P1[2] * (- one) * ( +cI * ( V3[2] + V3[5] ) ) + P1[3] * ( -V3[3] + cI * ( V3[4] ) ) ) ) ) + ( F2[5] * ( P1[0] * ( -V3[2] + V3[5] ) + ( P1[1] * ( V3[3] + cI * ( V3[4] ) ) + ( P1[2] * ( - cI * ( V3[3] ) + V3[4] ) + P1[3] * ( -V3[2] + V3[5] ) ) ) ) + M1 * ( F2[2] * ( -V3[3] + cI * ( V3[4] ) ) + F2[3] * ( V3[2] + V3[5] ) ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  // Compute the output wavefunction 'F2[6]' from the input wavefunctions F1[6], V3[6]
  __device__
  void FFV1_2( const cxtype_sv F1[],
               const cxtype_sv V3[],
               const cxtype COUP,
               const fptype M2,
               const fptype W2,
               cxtype_sv F2[] )
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype cI = cxmake( 0., 1. );
    F2[0] = + F1[0] + V3[0];
    F2[1] = + F1[1] + V3[1];
    const fptype_sv P2[4] = { -cxreal( F2[0] ), -cxreal( F2[1] ), -cximag( F2[1] ), -cximag( F2[0] ) };
    constexpr fptype one( 1. );
    const cxtype_sv denom = COUP / ( (P2[0] * P2[0] ) - ( P2[1] * P2[1] ) - ( P2[2] * P2[2] ) - ( P2[3] * P2[3] ) - M2 * ( M2 - cI * W2 ) );
    F2[2] = denom * cI * ( F1[2] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P2[2] * ( +cI * ( V3[3] ) - V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + ( F1[3] * ( P2[0] * ( V3[3]- cI * ( V3[4] ) ) + ( P2[1] * ( -V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] )- cI * ( V3[5] ) ) + P2[3] * ( -V3[3] + cI * ( V3[4] ) ) ) ) ) + M2 * ( F1[4] * ( V3[2] - V3[5] ) + F1[5] * ( -V3[3] + cI * ( V3[4] ) ) ) ) );
    F2[3] = denom * (- cI) * ( F1[2] * ( P2[0] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P2[1] * ( V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] + V3[5] ) ) - P2[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + ( F1[3] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3]- cI * ( V3[4] ) ) + ( P2[2] * ( +cI * ( V3[3] ) + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + M2 * ( F1[4] * ( V3[3] + cI * ( V3[4] ) ) - F1[5] * ( V3[2] + V3[5] ) ) ) );
    F2[4] = denom * (- cI) * ( F1[4] * ( P2[0] * ( -V3[2] + V3[5] ) + ( P2[1] * ( V3[3] + cI * ( V3[4] ) ) + ( P2[2] * ( - cI * ( V3[3] ) + V3[4] ) + P2[3] * ( -V3[2] + V3[5] ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[3]- cI * ( V3[4] ) ) + ( P2[1] * (- one) * ( V3[2] + V3[5] ) + ( P2[2] * ( +cI * ( V3[2] + V3[5] ) ) + P2[3] * ( V3[3]- cI * ( V3[4] ) ) ) ) ) + M2 * ( F1[2] * (- one) * ( V3[2] + V3[5] ) + F1[3] * ( -V3[3] + cI * ( V3[4] ) ) ) ) );
    F2[5] = denom * cI * ( F1[4] * ( P2[0] * (- one) * ( V3[3] + cI * ( V3[4] ) ) + ( P2[1] * ( V3[2] - V3[5] ) + ( P2[2] * ( +cI * ( V3[2] )- cI * ( V3[5] ) ) + P2[3] * ( V3[3] + cI * ( V3[4] ) ) ) ) ) + ( F1[5] * ( P2[0] * ( V3[2] + V3[5] ) + ( P2[1] * ( -V3[3] + cI * ( V3[4] ) ) + ( P2[2] * (- one) * ( +cI * ( V3[3] ) + V3[4] ) - P2[3] * ( V3[2] + V3[5] ) ) ) ) + M2 * ( F1[2] * ( V3[3] + cI * ( V3[4] ) ) + F1[3] * ( V3[2] - V3[5] ) ) ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

} // end namespace

