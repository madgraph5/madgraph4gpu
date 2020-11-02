//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

namespace MG5_sm
{
  using mgOnGpu::nw6;

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC inline const fptype& pIparIp4Ievt( T_Acc const &acc, const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    using mgOnGpu::np4;
#if defined MGONGPU_LAYOUT_ASA
    using mgOnGpu::npar;
    using mgOnGpu::nepp;
    const int ipag = ievt/nepp; // #eventpage in this iteration
    const int iepp = ievt%nepp; // #event in the current eventpage in this iteration
    // ASA: allmomenta[npag][npar][np4][nepp]
    return allmomenta[ipag*npar*np4*nepp + ipar*nepp*np4 + ip4*nepp + iepp]; // AOSOA[ipag][ipar][ip4][iepp]
#elif defined MGONGPU_LAYOUT_SOA
    const int nevt = blockDim.x * gridDim.x;
    // SOA: allmomenta[npar][np4][ndim]
    return allmomenta[ipar*np4*nevt + ip4*nevt + ievt]; // SOA[ipar][ip4][ievt]
#elif defined MGONGPU_LAYOUT_AOS
    using mgOnGpu::npar;
    // AOS: allmomenta[ndim][npar][np4]
    return allmomenta[ievt*npar*np4 + ipar*np4 + ip4]; // AOS[ievt][ipar][ip4]
#endif
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void imzxxxM0( T_Acc const &acc, const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
#if defined MGONGPU_WFMEM_LOCAL
                 cxtypeparam fisp[nw6],
#else
                 cxtypeparam fivp[],             // output: wavefunctions[5 * 6 * #threads_in_block]
#endif
                 const int ipar )          // input: particle# out of npar
  {
    {
#if !defined MGONGPU_WFMEM_LOCAL
      const int neib = blockDim.x; // number of events (threads) in block
      const int ieib = threadIdx.x; // index of event (thread) in block
#endif
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "imzxxxM0: ievt=%d ieib=%d\n", ievt, threadIdx.x );

      const fptype& pvec0 = pIparIp4Ievt( acc, allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( acc, allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( acc, allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( acc, allmomenta, ipar, 3, ievt );
#if !defined MGONGPU_WFMEM_LOCAL
      cxtyperef fi_0(acc, fivp[ipar*nw6*neib + 0*neib + ieib]);
      cxtyperef fi_1(acc, fivp[ipar*nw6*neib + 1*neib + ieib]);
      cxtyperef fi_2(acc, fivp[ipar*nw6*neib + 2*neib + ieib]);
      cxtyperef fi_3(acc, fivp[ipar*nw6*neib + 3*neib + ieib]);
      cxtyperef fi_4(acc, fivp[ipar*nw6*neib + 4*neib + ieib]);
      cxtyperef fi_5(acc, fivp[ipar*nw6*neib + 5*neib + ieib]);
#else
      cxtyperef fi_0(acc, fisp[0]);
      cxtyperef fi_1(acc, fisp[1]);
      cxtyperef fi_2(acc, fisp[2]);
      cxtyperef fi_3(acc, fisp[3]);
      cxtyperef fi_4(acc, fisp[4]);
      cxtyperef fi_5(acc, fisp[5]);
#endif
      fi_0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi_1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX = PY = 0 and E = -P3 > 0)
      {
        cxtypeparam chip[2] = { cxmake( 0, 0 ), cxmake( -nhel * sqrt(acc, 2 * pvec0), 0 ) };
        cxtyperef chi0(acc, chip[0]);
        cxtyperef chi1(acc, chip[1]);
        if (nh == 1)
        {
          fi_2 = cxmake( 0, 0 );
          fi_3 = cxmake( 0, 0 );
          fi_4 = chi0;
          fi_5 = chi1;
        }
        else
        {
          fi_2 = chi1;
          fi_3 = chi0;
          fi_4 = cxmake( 0, 0 );
          fi_5 = cxmake( 0, 0 );
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void ixzxxxM0( T_Acc const &acc, const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
#if defined MGONGPU_WFMEM_LOCAL
                 cxtypeparam fisp[nw6],
#else
                 cxtypeparam fivp[],             // output: wavefunctions[5 * 6 * #threads_in_block]
#endif
                 const int ipar )          // input: particle# out of npar
  {
    {
#if !defined MGONGPU_WFMEM_LOCAL
      const int neib = blockDim.x; // number of events (threads) in block
      const int ieib = threadIdx.x; // index of event (thread) in block
#endif
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixzxxxM0: ievt=%d ieib=%d\n", ievt, threadIdx.x );

      const fptype& pvec0 = pIparIp4Ievt( acc, allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( acc, allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( acc, allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( acc, allmomenta, ipar, 3, ievt );
#if !defined MGONGPU_WFMEM_LOCAL
      cxtyperef fi_0(acc, fivp[ipar*nw6*neib + 0*neib + ieib]);
      cxtyperef fi_1(acc, fivp[ipar*nw6*neib + 1*neib + ieib]);
      cxtyperef fi_2(acc, fivp[ipar*nw6*neib + 2*neib + ieib]);
      cxtyperef fi_3(acc, fivp[ipar*nw6*neib + 3*neib + ieib]);
      cxtyperef fi_4(acc, fivp[ipar*nw6*neib + 4*neib + ieib]);
      cxtyperef fi_5(acc, fivp[ipar*nw6*neib + 5*neib + ieib]);
#else
      cxtyperef fi_0(acc, fisp[0]);
      cxtyperef fi_1(acc, fisp[1]);
      cxtyperef fi_2(acc, fisp[2]);
      cxtyperef fi_3(acc, fisp[3]);
      cxtyperef fi_4(acc, fisp[4]);
      cxtyperef fi_5(acc, fisp[5]);
#endif
      fi_0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi_1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX and PY are not 0)
      {
        const fptype sqp0p3 = sqrt( acc, pvec0 + pvec3 ) * nsf;
        cxtypeparam chip[2] = { cxmake( sqp0p3, 0 ), cxmake( nh * pvec1 / sqp0p3, pvec2 / sqp0p3 ) };
        cxtyperef chi0(acc, chip[0]);
        cxtyperef chi1(acc, chip[1]);
        if ( nh == 1 )
        {
          fi_2 = cxmake( 0, 0 );
          fi_3 = cxmake( 0, 0 );
          fi_4 = chi0;
          fi_5 = chi1;
        }
        else
        {
          fi_2 = chi1;
          fi_3 = chi0;
          fi_4 = cxmake( 0, 0 );
          fi_5 = cxmake( 0, 0 );
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void oxzxxxM0( T_Acc const &acc, const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
#if defined MGONGPU_WFMEM_LOCAL
                 cxtypeparam fosp[nw6],
#else
                 cxtypeparam fovp[],             // output: wavefunctions[5 * 6 * #threads_in_block]
#endif
                 const int ipar )          // input: particle# out of npar
  {
    {
#if !defined MGONGPU_WFMEM_LOCAL
      const int neib = blockDim.x; // number of events (threads) in block
      const int ieib = threadIdx.x; // index of event (thread) in block
#endif
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxzxxxM0: ievt=%d ieib=%d\n", ievt, threadIdx.x );

      const fptype& pvec0 = pIparIp4Ievt( acc, allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( acc, allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( acc, allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( acc, allmomenta, ipar, 3, ievt );
// dhs
//  std::cout << "oxzxxxM0 pvec[0-3]=" << pvec0 << " " << pvec1 << " " << pvec2 << " " << pvec3 << std::endl;

#if !defined MGONGPU_WFMEM_LOCAL
      cxtyperef fo_0(acc, fovp[ipar*nw6*neib + 0*neib + ieib]);
      cxtyperef fo_1(acc, fovp[ipar*nw6*neib + 1*neib + ieib]);
      cxtyperef fo_2(acc, fovp[ipar*nw6*neib + 2*neib + ieib]);
      cxtyperef fo_3(acc, fovp[ipar*nw6*neib + 3*neib + ieib]);
      cxtyperef fo_4(acc, fovp[ipar*nw6*neib + 4*neib + ieib]);
      cxtyperef fo_5(acc, fovp[ipar*nw6*neib + 5*neib + ieib]);
#else
      cxtyperef fo_0(acc, fosp[0]);
      cxtyperef fo_1(acc, fosp[1]);
      cxtyperef fo_2(acc, fosp[2]);
      cxtyperef fo_3(acc, fosp[3]);
      cxtyperef fo_4(acc, fosp[4]);
      cxtyperef fo_5(acc, fosp[5]);
#endif
      fo_0 = cxmake( pvec0 * nsf, pvec3 * nsf );
      fo_1 = cxmake( pvec1 * nsf, pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // EITHER (Px and Py are not zero)
      // OR (PX = PY = 0 and E = P3 > 0)
      {
        const fptype sqp0p3 = sqrt( acc, pvec0 + pvec3 ) * nsf;
        cxtypeparam chip[2] = { cxmake( sqp0p3, 0 ), cxmake( nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 ) };
        cxtyperef chi0(acc, chip[0]);
        cxtyperef chi1(acc, chip[1]);
        if( nh == 1 )
        {
          fo_2 = chi0;
          fo_3 = chi1;
          fo_4 = cxmake( 0, 0 );
          fo_5 = cxmake( 0, 0 );
        }
        else
        {
          fo_2 = cxmake( 0, 0 );
          fo_3 = cxmake( 0, 0 );
          fo_4 = chi1;
          fo_5 = chi0;
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void FFV1_0( T_Acc const &acc,
#if !defined MGONGPU_WFMEM_LOCAL
              const cxtypeparam F1Vp[],  // input: wavefunctionA[6 * #threads_in_block]
              const cxtypeparam F2Vp[],  // input: wavefunctionB[6 * #threads_in_block]
              const cxtypeparam V3Vp[],  // input: wavefunctionC[6 * #threads_in_block]
#else
              const cxtypeparam F1Sp[],   // input wavefunctionA[6]
              const cxtypeparam F2Sp[],   // input wavefunctionB[6]
              const cxtypeparam V3Sp[],   // input wavefunctionC[6]
#endif
              const cxtypeparam COUPpr,
              cxtypeparam* vertexp )    // output
  {
    const cxtypeconstref COUP(acc, COUPpr);
    cxtyperef vertex(acc, *vertexp);

#if !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtypeconstref F1_2(acc, F1Vp[2*neib + ieib]);
    const cxtypeconstref F1_3(acc, F1Vp[3*neib + ieib]);
    const cxtypeconstref F1_4(acc, F1Vp[4*neib + ieib]);
    const cxtypeconstref F1_5(acc, F1Vp[5*neib + ieib]);
    const cxtypeconstref F2_2(acc, F2Vp[2*neib + ieib]);
    const cxtypeconstref F2_3(acc, F2Vp[3*neib + ieib]);
    const cxtypeconstref F2_4(acc, F2Vp[4*neib + ieib]);
    const cxtypeconstref F2_5(acc, F2Vp[5*neib + ieib]);
    const cxtypeconstref V3_2(acc, V3Vp[2*neib + ieib]);
    const cxtypeconstref V3_3(acc, V3Vp[3*neib + ieib]);
    const cxtypeconstref V3_4(acc, V3Vp[4*neib + ieib]);
    const cxtypeconstref V3_5(acc, V3Vp[5*neib + ieib]);
#else
    const cxtypeconstref F1_2(acc, F1Sp[2]);
    const cxtypeconstref F1_3(acc, F1Sp[3]);
    const cxtypeconstref F1_4(acc, F1Sp[4]);
    const cxtypeconstref F1_5(acc, F1Sp[5]);
    const cxtypeconstref F2_2(acc, F2Sp[2]);
    const cxtypeconstref F2_3(acc, F2Sp[3]);
    const cxtypeconstref F2_4(acc, F2Sp[4]);
    const cxtypeconstref F2_5(acc, F2Sp[5]);
    const cxtypeconstref V3_2(acc, V3Sp[2]);
    const cxtypeconstref V3_3(acc, V3Sp[3]);
    const cxtypeconstref V3_4(acc, V3Sp[4]);
    const cxtypeconstref V3_5(acc, V3Sp[5]);
#endif
    cxtypeparam cxp[2];
    cxtyperef cI(acc, cxp[0]);
    cxtyperef TMP4(acc, cxp[1]);
    cI = cxmake( 0, 1 );
    TMP4 =
      ( F1_2 * ( F2_4 * ( V3_2 + V3_5 ) +
                 F2_5 * ( V3_3 + cI * ( V3_4 ) ) ) +
        ( F1_3 * ( F2_4 * ( V3_3 - cI * ( V3_4) )
                   + F2_5 * ( V3_2 - V3_5 ) ) +
          ( F1_4 * ( F2_2 * ( V3_2 - V3_5 )
                     - F2_3 * ( V3_3 + cI * ( V3_4 ) ) ) +
            F1_5 * ( F2_2 * ( -V3_3 + cI * ( V3_4 ) )
                     + F2_3 * ( V3_2 + V3_5 ) )
            )
          )
        );
    vertex = COUP * ( -cI ) * TMP4;
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void FFV1P0_3( T_Acc const &acc,
#if !defined MGONGPU_WFMEM_LOCAL
                const cxtypeparam F1Vp[],  // input: wavefunctionA[6 * #threads_in_block]
                const cxtypeparam F2Vp[],  // input: wavefunctionB[6 * #threads_in_block]
#else
                const cxtypeparam F1Sp[],   // input wavefunctionA[6]
                const cxtypeparam F2Sp[],   // input wavefunctionB[6]
#endif
                const cxtypeparam COUPpr,
                const fptype M3,
                const fptype W3
#if !defined MGONGPU_WFMEM_LOCAL
                , cxtypeparam V3Vp[]         // output: wavefunctionC[6 * #threads_in_block]
#else
                , cxtypeparam V3Sp[]          // output wavefunctionC[6]
#endif
                )
  {
    const cxtypeconstref COUP(acc, COUPpr);

#if !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtypeconstref F1_0(acc, F1Vp[0*neib + ieib]);
    const cxtypeconstref F1_1(acc, F1Vp[1*neib + ieib]);
    const cxtypeconstref F1_2(acc, F1Vp[2*neib + ieib]);
    const cxtypeconstref F1_3(acc, F1Vp[3*neib + ieib]);
    const cxtypeconstref F1_4(acc, F1Vp[4*neib + ieib]);
    const cxtypeconstref F1_5(acc, F1Vp[5*neib + ieib]);
    const cxtypeconstref F2_0(acc, F2Vp[0*neib + ieib]);
    const cxtypeconstref F2_1(acc, F2Vp[1*neib + ieib]);
    const cxtypeconstref F2_2(acc, F2Vp[2*neib + ieib]);
    const cxtypeconstref F2_3(acc, F2Vp[3*neib + ieib]);
    const cxtypeconstref F2_4(acc, F2Vp[4*neib + ieib]);
    const cxtypeconstref F2_5(acc, F2Vp[5*neib + ieib]);
    cxtyperef V3_0(acc, V3Vp[0*neib + ieib]);
    cxtyperef V3_1(acc, V3Vp[1*neib + ieib]);
    cxtyperef V3_2(acc, V3Vp[2*neib + ieib]);
    cxtyperef V3_3(acc, V3Vp[3*neib + ieib]);
    cxtyperef V3_4(acc, V3Vp[4*neib + ieib]);
    cxtyperef V3_5(acc, V3Vp[5*neib + ieib]);
#else
    const cxtypeconstref F1_0(acc, F1Sp[0]);
    const cxtypeconstref F1_1(acc, F1Sp[1]);
    const cxtypeconstref F1_2(acc, F1Sp[2]);
    const cxtypeconstref F1_3(acc, F1Sp[3]);
    const cxtypeconstref F1_4(acc, F1Sp[4]);
    const cxtypeconstref F1_5(acc, F1Sp[5]);
    const cxtypeconstref F2_0(acc, F2Sp[0]);
    const cxtypeconstref F2_1(acc, F2Sp[1]);
    const cxtypeconstref F2_2(acc, F2Sp[2]);
    const cxtypeconstref F2_3(acc, F2Sp[3]);
    const cxtypeconstref F2_4(acc, F2Sp[4]);
    const cxtypeconstref F2_5(acc, F2Sp[5]);
    cxtyperef V3_0(acc, V3Sp[0]);
    cxtyperef V3_1(acc, V3Sp[1]);
    cxtyperef V3_2(acc, V3Sp[2]);
    cxtyperef V3_3(acc, V3Sp[3]);
    cxtyperef V3_4(acc, V3Sp[4]);
    cxtyperef V3_5(acc, V3Sp[5]);
#endif
    cxtypeparam cxp[2];
    cxtyperef cI(acc, cxp[0]);
    cxtyperef denom(acc, cxp[1]);
    cI = cxmake( 0, 1 );
    V3_0 = + F1_0 + F2_0;
    V3_1 = + F1_1 + F2_1;
    const fptype PPP0 = -cxreal( V3_0 );
    const fptype PPP1 = -cxreal( V3_1 );
    const fptype PPP2 = -cximag( V3_1 );
    const fptype PPP3 = -cximag( V3_0 );
//dhs
//    denom = COUP;
//std::cout << "FFV1P0_3 denom_top=(" << cxp[1].r << "," << cxp[1].i << ")" << std::endl;
//    denom = ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) - ( PPP3 * PPP3 ) - M3 * ( M3 - cI * W3 ) );
//std::cout << "FFV1P0_3 denom_bot=(" << cxp[1].r << "," << cxp[1].i << ")" << std::endl;
    denom =
      COUP / ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) - ( PPP3 * PPP3 ) - M3 * ( M3 - cI * W3 ) );
//std::cout << "Done div" << std::endl;
    V3_2 = denom * ( -cI ) * ( F1_2 * F2_4 + F1_3 * F2_5 + F1_4 * F2_2 + F1_5 * F2_3 );
//dhs
//std::cout << "FFV1P0_3 denom=(" << cxp[1].r << "," << cxp[1].i << ")" << std::endl;
    V3_3 = denom * ( -cI ) * ( -F1_2 * F2_5 - F1_3 * F2_4 + F1_4 * F2_3 + F1_5 * F2_2 );
//dhs
//std::cout << "FFV1P0_3 V3_3=(" << V3Sp[3].r << "," << V3Sp[3].i << ")" << std::endl;
    V3_4 = denom * ( -cI ) * ( -cI * ( F1_2 * F2_5 + F1_5 * F2_2 ) + cI * ( F1_3 * F2_4 + F1_4 * F2_3 ) );
    V3_5 = denom * ( -cI ) * ( -F1_2 * F2_4 - F1_5 * F2_3 + F1_3 * F2_5 + F1_4 * F2_2 );
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void FFV2_4_0( T_Acc const &acc,
#if !defined MGONGPU_WFMEM_LOCAL
                const cxtypeparam F1Vp[],  // input: wavefunctionA[6 * #threads_in_block]
                const cxtypeparam F2Vp[],  // input: wavefunctionB[6 * #threads_in_block]
                const cxtypeparam V3Vp[],  // input: wavefunctionC[6 * #threads_in_block]
#else
                const cxtypeparam F1Sp[],   // input wavefunctionA[6]
                const cxtypeparam F2Sp[],   // input wavefunctionB[6]
                const cxtypeparam V3Sp[],   // input wavefunctionC[6]
#endif
                const cxtypeparam COUP1pr,
                const cxtypeparam COUP2pr,
                cxtypeparam* vertexp )    // output
  {
    const cxtypeconstref COUP1(acc, COUP1pr);
    const cxtypeconstref COUP2(acc, COUP2pr);
    cxtyperef vertex(acc, *vertexp);

#if !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtypeconstref F1_2(acc, F1Vp[2*neib + ieib]);
    const cxtypeconstref F1_3(acc, F1Vp[3*neib + ieib]);
    const cxtypeconstref F1_4(acc, F1Vp[4*neib + ieib]);
    const cxtypeconstref F1_5(acc, F1Vp[5*neib + ieib]);
    const cxtypeconstref F2_2(acc, F2Vp[2*neib + ieib]);
    const cxtypeconstref F2_3(acc, F2Vp[3*neib + ieib]);
    const cxtypeconstref F2_4(acc, F2Vp[4*neib + ieib]);
    const cxtypeconstref F2_5(acc, F2Vp[5*neib + ieib]);
    const cxtypeconstref V3_2(acc, V3Vp[2*neib + ieib]);
    const cxtypeconstref V3_3(acc, V3Vp[3*neib + ieib]);
    const cxtypeconstref V3_4(acc, V3Vp[4*neib + ieib]);
    const cxtypeconstref V3_5(acc, V3Vp[5*neib + ieib]);
#else
    const cxtypeconstref F1_2(acc, F1Sp[2]);
    const cxtypeconstref F1_3(acc, F1Sp[3]);
    const cxtypeconstref F1_4(acc, F1Sp[4]);
    const cxtypeconstref F1_5(acc, F1Sp[5]);
    const cxtypeconstref F2_2(acc, F2Sp[2]);
    const cxtypeconstref F2_3(acc, F2Sp[3]);
    const cxtypeconstref F2_4(acc, F2Sp[4]);
    const cxtypeconstref F2_5(acc, F2Sp[5]);
    const cxtypeconstref V3_2(acc, V3Sp[2]);
    const cxtypeconstref V3_3(acc, V3Sp[3]);
    const cxtypeconstref V3_4(acc, V3Sp[4]);
    const cxtypeconstref V3_5(acc, V3Sp[5]);
#endif
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    cxtypeparam cxp[3];
    cxtyperef cI(acc, cxp[0]);
    cxtyperef TMP2(acc, cxp[1]);
    cxtyperef TMP0(acc, cxp[2]);
    cI = cxmake( 0, 1 );
    TMP2 =
      ( F1_4 * ( F2_2 * ( V3_2 - V3_5 ) - F2_3 * ( V3_3 + cI * ( V3_4 ) ) ) +
        F1_5 * ( F2_2 * ( -V3_3 + cI * ( V3_4 ) ) + F2_3 * ( V3_2 + V3_5 ) ) );
    TMP0 =
      ( F1_2 * ( F2_4 * ( V3_2 + V3_5 ) + F2_5 * ( V3_3 + cI * ( V3_4 ) ) ) +
        F1_3 * ( F2_4 * ( V3_3 - cI * ( V3_4 ) ) + F2_5 * ( V3_2 - V3_5 ) ) );
    vertex = -fp1 * ( COUP2 * ( +cI * ( TMP0 ) + fp2 * cI * ( TMP2 ) ) + cI * ( TMP0 * COUP1 ) );
  }

  //--------------------------------------------------------------------------

  template< typename T_Acc >
  ALPAKA_FN_ACC void FFV2_4_3( T_Acc const &acc,
#if !defined MGONGPU_WFMEM_LOCAL
                const cxtypeparam F1Vp[],  // input: wavefunctionA[6 * #threads_in_block]
                const cxtypeparam F2Vp[],  // input: wavefunctionB[6 * #threads_in_block]
#else
                const cxtypeparam F1Sp[],   // input wavefunctionA[6]
                const cxtypeparam F2Sp[],   // input wavefunctionB[6]
#endif
                const cxtypeparam COUP1pr,
                const cxtypeparam COUP2pr,
                const fptype M3,
                const fptype W3
#if !defined MGONGPU_WFMEM_LOCAL
                , cxtypeparam V3Vp[]         // output: wavefunctionC[6 * #threads_in_block]
#else
                , cxtypeparam V3Sp[]          // output wavefunctionC[6]
#endif
                )
  {
    const cxtypeconstref COUP1(acc, COUP1pr);
    const cxtypeconstref COUP2(acc, COUP2pr);

#if !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtypeconstref F1_0(acc, F1Vp[0*neib + ieib]);
    const cxtypeconstref F1_1(acc, F1Vp[1*neib + ieib]);
    const cxtypeconstref F1_2(acc, F1Vp[2*neib + ieib]);
    const cxtypeconstref F1_3(acc, F1Vp[3*neib + ieib]);
    const cxtypeconstref F1_4(acc, F1Vp[4*neib + ieib]);
    const cxtypeconstref F1_5(acc, F1Vp[5*neib + ieib]);
    const cxtypeconstref F2_0(acc, F2Vp[0*neib + ieib]);
    const cxtypeconstref F2_1(acc, F2Vp[1*neib + ieib]);
    const cxtypeconstref F2_2(acc, F2Vp[2*neib + ieib]);
    const cxtypeconstref F2_3(acc, F2Vp[3*neib + ieib]);
    const cxtypeconstref F2_4(acc, F2Vp[4*neib + ieib]);
    const cxtypeconstref F2_5(acc, F2Vp[5*neib + ieib]);
    cxtyperef V3_0(acc, V3Vp[0*neib + ieib]);
    cxtyperef V3_1(acc, V3Vp[1*neib + ieib]);
    cxtyperef V3_2(acc, V3Vp[2*neib + ieib]);
    cxtyperef V3_3(acc, V3Vp[3*neib + ieib]);
    cxtyperef V3_4(acc, V3Vp[4*neib + ieib]);
    cxtyperef V3_5(acc, V3Vp[5*neib + ieib]);
#else
    const cxtypeconstref F1_0(acc, F1Sp[0]);
    const cxtypeconstref F1_1(acc, F1Sp[1]);
    const cxtypeconstref F1_2(acc, F1Sp[2]);
    const cxtypeconstref F1_3(acc, F1Sp[3]);
    const cxtypeconstref F1_4(acc, F1Sp[4]);
    const cxtypeconstref F1_5(acc, F1Sp[5]);
    const cxtypeconstref F2_0(acc, F2Sp[0]);
    const cxtypeconstref F2_1(acc, F2Sp[1]);
    const cxtypeconstref F2_2(acc, F2Sp[2]);
    const cxtypeconstref F2_3(acc, F2Sp[3]);
    const cxtypeconstref F2_4(acc, F2Sp[4]);
    const cxtypeconstref F2_5(acc, F2Sp[5]);
    cxtyperef V3_0(acc, V3Sp[0]);
    cxtyperef V3_1(acc, V3Sp[1]);
    cxtyperef V3_2(acc, V3Sp[2]);
    cxtyperef V3_3(acc, V3Sp[3]);
    cxtyperef V3_4(acc, V3Sp[4]);
    cxtyperef V3_5(acc, V3Sp[5]);
#endif
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    cxtypeparam cxp[4];
    cxtyperef cI(acc, cxp[0]);
    cxtyperef TMP1(acc, cxp[1]);
    cxtyperef TMP3(acc, cxp[2]);
    cxtyperef denom(acc, cxp[3]);
    cI = cxmake( 0, 1 );
    fptype OM3 = 0;
    if ( M3 != 0 ) OM3 = fp1 / ( M3 * M3 );
    V3_0 = + F1_0 + F2_0;
    V3_1 = + F1_1 + F2_1;
    const fptype PPP0 = -cxreal( V3_0 );
    const fptype PPP1 = -cxreal( V3_1 );
    const fptype PPP2 = -cximag( V3_1 );
    const fptype PPP3 = -cximag( V3_0 );
    TMP1 =
      ( F1_2 * ( F2_4 * ( PPP0 + PPP3 ) + F2_5 * ( PPP1 + cI * ( PPP2 ) ) ) +
        F1_3 * ( F2_4 * ( PPP1 - cI * ( PPP2 ) ) + F2_5 * ( PPP0 - PPP3 ) ) );
    TMP3 =
      ( F1_4 * ( F2_2 * ( PPP0 - PPP3 ) - F2_3 * ( PPP1 + cI * ( PPP2 ) ) ) +
        F1_5 * ( F2_2 * ( -PPP1 + cI * ( PPP2 ) ) + F2_3 * ( PPP0 + PPP3 ) ) );
    denom =
      fp1 / ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) -
              ( PPP3 * PPP3 ) - M3 * (M3 - cI * W3 ) );
    V3_2 = denom * ( -fp2 * cI ) *
      ( COUP2 * ( OM3 * ( -fp1/fp2 ) * PPP0 * ( TMP1 + fp2 * ( TMP3 ) )
                  + ( +fp1/fp2 * ( F1_2 * F2_4 + F1_3 * F2_5 ) + F1_4 * F2_2 + F1_5 * F2_3 ) )
        + fp1/fp2 * ( COUP1 * ( F1_2 * F2_4 + F1_3 * F2_5 - PPP0 * OM3 * TMP1 ) ) );
    V3_3 = denom * ( -fp2 * cI ) *
      ( COUP2 * ( OM3 * ( -fp1/fp2 ) * PPP1 * ( TMP1 + fp2 * ( TMP3 ) )
                  + ( -fp1/fp2 * ( F1_2 * F2_5 + F1_3 * F2_4 ) + F1_4 * F2_3 + F1_5 * F2_2 ) )
        - fp1/fp2 * ( COUP1 * ( F1_2 * F2_5 + F1_3 * F2_4 + PPP1 * OM3 * TMP1 ) ) );
    V3_4 = denom * cI *
      ( COUP2 * ( OM3 * PPP2 * ( TMP1 + fp2 * ( TMP3 ) )
                  + ( +cI * ( F1_2 * F2_5 ) - cI * ( F1_3 * F2_4 )
                      - fp2 * cI * ( F1_4 * F2_3 )
                      + fp2 * cI * ( F1_5 * F2_2 ) ) )
        + COUP1 * ( +cI * ( F1_2 * F2_5 ) - cI * ( F1_3 * F2_4 ) + PPP2 * OM3 * TMP1 ) );
    V3_5 = denom * fp2 * cI *
      ( COUP2 * ( OM3 * fp1/fp2 * PPP3 * ( TMP1 + fp2 * ( TMP3 ) ) +
                  ( +fp1/fp2 * ( F1_2 * F2_4 ) - fp1/fp2 * ( F1_3 * F2_5 ) - F1_4 * F2_2 + F1_5 * F2_3 ) )
        + fp1/fp2 * ( COUP1 * ( F1_2 * F2_4 + PPP3 * OM3 * TMP1 - F1_3 * F2_5 ) ) );
  }


}  // end namespace


//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <algorithm>
#include <iostream>

#include "mgOnGpuConfig.h"

#include "CPPProcess.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

namespace gProc
{
  using mgOnGpu::np4;
  using mgOnGpu::npar;
  const int ncomb = 16; // #helicity combinations is hardcoded for this process (eemumu: ncomb=16)

  ALPAKA_STATIC_ACC_MEM_CONSTANT int cHel[ncomb][npar];
  ALPAKA_STATIC_ACC_MEM_CONSTANT fptype cIPC[6];  // coupling ?
  ALPAKA_STATIC_ACC_MEM_CONSTANT fptype cIPD[2];

  ALPAKA_STATIC_ACC_MEM_GLOBAL unsigned long long sigmakin_itry = 0; // first iteration over nevt events
  ALPAKA_STATIC_ACC_MEM_GLOBAL bool sigmakin_goodhel[ncomb] = { false };

  //--------------------------------------------------------------------------

  using mgOnGpu::nwf;
  using mgOnGpu::nw6;

#if defined MGONGPU_WFMEM_SHARED

  int sigmakin_sharedmem_nbytes( const int ntpb ) // input: #threads per block
  {
    // Wavefunctions for this block: cxtype bwf[5 * 6 * #threads_in_block]
    const int nbytesBwf = nwf * nw6 * ntpb * sizeof(cxtypeparam);
    return nbytesBwf;
  }

#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess

  // ** NB: allmomenta can have three different layouts
  // ASA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
  // SOA: allmomenta[npar][np4][ndim]
  // AOS: allmomenta[ndim][npar][np4]

  template< typename T_Acc >
  ALPAKA_FN_ACC void calculate_wavefunctions( T_Acc const &acc, int ihel,
                                const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
    fptype &matrix
#if defined MGONGPU_WFMEM_GLOBAL
                                , cxtypeparam* tmpWFsp          // tmp[(nwf=5)*(nw6=6)*(nevt=nblk*ntpb)] 
#endif
                                )
  {
    cxtypeparam ampp[2];
#if !defined MGONGPU_WFMEM_LOCAL
    // eventually move to same AOSOA everywhere, blocks and threads
#if defined MGONGPU_WFMEM_GLOBAL
    const int iblk = blockIdx.x; // index of block in grid
    const int ntpb = blockDim.x; // index of block in grid
    cxtypeparam* bwf = &tmpWFsp[iblk*nwf*nw6*ntpb];
#elif defined MGONGPU_WFMEM_SHARED
    sharedMemExtern(bwf, cxtypeparam);
#endif
    MG5_sm::oxzxxxM0( acc, allmomenta, cHel[ihel][0], -1, bwf, 0 );
    MG5_sm::imzxxxM0( acc, allmomenta, cHel[ihel][1], +1, bwf, 1 );
    MG5_sm::ixzxxxM0( acc, allmomenta, cHel[ihel][2], -1, bwf, 2 );
    MG5_sm::oxzxxxM0( acc, allmomenta, cHel[ihel][3], +1, bwf, 3 );
#else //local
    cxtypeparam w[nwf][nw6]; // w[5][6]

// dhs
//#include <mutex>
//static std::mutex mdhs;
//std::lock_guard<std::mutex> mgu(mdhs);

/*
  std::cout << "<<< (blockIdx.x=" << blockIdx.x << ",threadIdx.x=" << threadIdx.x << " ihel=" << ihel << ") w in: " << std::endl;
  for (int i = 0; i < 5; ++i) {
    std::cout << "w" << i << ": ";
    for (int j = 0; j < 6; ++j) {
      if (w[i][j].r || w[i][j].i)
        std::cout << "(" << w[i][j].r << "," << w[i][j].i << ") ";
      else
        std::cout << "0 ";
    }
    std::cout << std::endl;
  }
*/

    MG5_sm::oxzxxxM0( acc, allmomenta, cHel[ihel][0], -1, w[0], 0 );
    MG5_sm::imzxxxM0( acc, allmomenta, cHel[ihel][1], +1, w[1], 1 );
    MG5_sm::ixzxxxM0( acc, allmomenta, cHel[ihel][2], -1, w[2], 2 );
    MG5_sm::oxzxxxM0( acc, allmomenta, cHel[ihel][3], +1, w[3], 3 );
#endif

#if !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block

    // Diagram 1
    MG5_sm::FFV1P0_3( acc,
                      &(bwf[1*nw6*neib]),
                      &(bwf[0*nw6*neib]),
                      cxmake( cIPC[0], cIPC[1] ),
                      0.,
                      0.,
                      &(bwf[4*nw6*neib]) );

    MG5_sm::FFV1_0( acc,
                    &(bwf[2*nw6*neib]),
                    &(bwf[3*nw6*neib]),
                    &(bwf[4*nw6*neib]),
                    cxmake( cIPC[0], cIPC[1] ),
                    &ampp[0] );

    // Diagram 2
    MG5_sm::FFV2_4_3( acc,
                      &(bwf[1*nw6*neib]),
                      &(bwf[0*nw6*neib]),
                      cxmake( cIPC[2], cIPC[3] ),
                      cxmake( cIPC[4], cIPC[5] ),
                      cIPD[0],
                      cIPD[1],
                      &(bwf[4*nw6*neib]) );

    MG5_sm::FFV2_4_0( acc,
                      &(bwf[2*nw6*neib]),
                      &(bwf[3*nw6*neib]),
                      &(bwf[4*nw6*neib]),
                      cxmake( cIPC[2], cIPC[3] ),
                      cxmake( cIPC[4], cIPC[5] ),
                      &ampp[1] );
#else

    // Diagram 1
    MG5_sm::FFV1P0_3( acc, w[1], w[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w[4] );
    MG5_sm::FFV1_0( acc, w[2], w[3], w[4], cxmake( cIPC[0], cIPC[1] ), &ampp[0] );
//std::cout << "using w[2],w[3],w[4] and (" << cIPC[0] << "," << cIPC[1] << ") got ampp[0]=(" << ampp[0].r << "," << ampp[0].i << ")" << std::endl;

    // Diagram 2
    MG5_sm::FFV2_4_3( acc, w[1], w[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w[4] );
    MG5_sm::FFV2_4_0( acc, w[2], w[3], w[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &ampp[1] );

// dhs
/*
  std::cout << "<<< w out: " << std::endl;
  for (int i = 0; i < 5; ++i) {
    std::cout << "w" << i << ": ";
    for (int j = 0; j < 6; ++j) {
      if (w[i][j].r || w[i][j].i)
        std::cout << "(" << w[i][j].r << "," << w[i][j].i << ") ";
      else
        std::cout << "0 ";
    }
    std::cout << std::endl;
    std::cout << "ampp[0]=(" << ampp[0].r << "," << ampp[0].i << ") (" << ampp[1].r << "," << ampp[1].i << ")";
    std::cout << std::endl;
  }
*/

#endif

    const int ncolor = 1;
    cxtypeparam jampp[ncolor];

    cxtypeptr jamp(acc, jampp);
    cxtypeptr amp(acc, ampp);

    // The color matrix;
    static const fptype denom[ncolor] = {1};
    static const fptype cf[ncolor][ncolor] = {{1}};

// dhs
//#include <mutex>
//static std::mutex mdhs;
//std::lock_guard<std::mutex> mgu(mdhs);

    // Calculate color flows
    jamp[0] = -amp[0] - amp[1];

    // Sum and square the color flows to get the matrix element
    for( int icol = 0; icol < ncolor; icol++ )
    {
      cxtypeparam zp = cxmake(0,0);
      cxtyperef ztemp(acc, zp);
      for( int jcol = 0; jcol < ncolor; jcol++ )
        ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
      matrix = matrix + cxreal( ztemp * conj( jamp[icol] ) ) / denom[icol];
// dhs
//std::cout << "threadIdx.x=" << threadIdx.x << " matrix=" << matrix << " jampp[0]=(" << jampp[0].r << "," << jampp[0].i << ") ampp[0]=(" << ampp[0].r << "," << ampp[0].i << ") ampp[1]=(" << ampp[1].r << "," << ampp[1].i << ")" << std::endl;
    }

    // Store the leading color flows for choice of color
    // for(i=0;i < ncolor; i++)
    // jamp2[0][i] += cxreal( jamp[i]*conj( jamp[i] ) );

  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( int numiterations,
                          int gpublocks,
                          int gputhreads,
                          bool verbose )
    : m_numiterations( numiterations )
    , gpu_nblocks( gpublocks )
    , gpu_nthreads( gputhreads )
    , dim( gpu_nblocks * gpu_nthreads )
    , m_verbose( verbose )
  {
    // Helicities for the process - nodim
    static const int tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };

    checkCupla( cuplaMemcpy( cHel, tHel, ncomb * nexternal * sizeof(int), cuplaMemcpyHostToDevice ) );

    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
//    assert( sizeof(cxtype) == 2*sizeof(fptype) );
  }

  //--------------------------------------------------------------------------

  CPPProcess::~CPPProcess() {}

  //--------------------------------------------------------------------------

  const std::vector<fptype> &CPPProcess::getMasses() const {return mME;}

  //--------------------------------------------------------------------------
  // Initialize process.

  void CPPProcess::initProc(std::string param_card_name)
  {
    // Instantiate the model class and set parameters that stay fixed during run
    pars = Parameters_sm::getInstance();
    SLHAReader slha(param_card_name, m_verbose);
    pars->setIndependentParameters(slha);
    pars->setIndependentCouplings();
    if (m_verbose) {
      pars->printIndependentParameters();
      pars->printIndependentCouplings();
    }
    pars->setDependentParameters();
    pars->setDependentCouplings();
    // Set external particle masses for this matrix element
    mME.push_back(pars->ZERO);
    mME.push_back(pars->ZERO);
    mME.push_back(pars->ZERO);
    mME.push_back(pars->ZERO);
    static cxtypeparam tIPC[3] = { cxmake( pars->GC_3 ), cxmake( pars->GC_50 ), cxmake( pars->GC_59 ) };
    static fptype tIPD[2] = { (fptype)pars->mdl_MZ, (fptype)pars->mdl_WZ };

    checkCupla( cuplaMemcpy( cIPC, tIPC, 3 * sizeof(cxtypeparam), cuplaMemcpyHostToDevice ) );
// dhs
//std::cout << "here00 tIPC[2]=(" << tIPC[2].r << "," << tIPC[2].i << ")" << std::endl;
    checkCupla( cuplaMemcpy( cIPD, tIPD, 2 * sizeof(fptype), cuplaMemcpyHostToDevice ) );

  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour.

  // ** NB: allmomenta can have three different layouts
  // ASA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
  // SOA: allmomenta[npar][np4][ndim]
  // AOS: allmomenta[ndim][npar][np4]
  template< typename T_Acc >
  ALPAKA_FN_ACC
  void sigmaKin::operator()( T_Acc const &acc,
               const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               fptype* output            // output[nevt]
               // NB: nevt == ndim=gpublocks*gputhreads
#if defined MGONGPU_WFMEM_GLOBAL
               , cxtypeparam* tmpWFsp          // tmp[(nwf=5)*(nw6=6)*nevt] 
#endif
               ) const
  {
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows
    const int maxtry = 10;

    {
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "sigmakin: ievt %d\n", ievt );

      // Denominators: spins, colors and identical particles
      const int nprocesses = 1;
      const int denominators[nprocesses] = {4};

      // Reset the matrix elements
      fptype matrix_element[nprocesses];
      for( int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] = 0.;
      }

      fptype melast = matrix_element[0];
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry>maxtry && !sigmakin_goodhel[ihel] ) continue;
        // Adds ME for ihel to matrix_element[0]

#if defined MGONGPU_WFMEM_GLOBAL
        calculate_wavefunctions(acc, ihel, allmomenta, matrix_element[0], tmpWFsp); 
#else
        calculate_wavefunctions(acc, ihel, allmomenta, matrix_element[0]); 
#endif
        if ( sigmakin_itry<=maxtry )
        {
          if ( !sigmakin_goodhel[ihel] && matrix_element[0]>melast ) sigmakin_goodhel[ihel] = true;
          melast = matrix_element[0];
        }
      }

      for (int iproc = 0; iproc < nprocesses; ++iproc)
      {
        matrix_element[iproc] /= denominators[iproc];
      }

      for (int iproc = 0; iproc < nprocesses; ++iproc)
      {
        output[iproc*nprocesses + ievt] = matrix_element[iproc];
      }

      if ( sigmakin_itry <= maxtry ) {
        atomicAdd(&sigmakin_itry, 1ull, ::alpaka::hierarchy::Grids{});
      }

    }
    // ** END LOOP ON IEVT **
  }

  //--------------------------------------------------------------------------

}
