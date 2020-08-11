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

#ifndef __CUDACC__
  // Quick and dirty way to share nevt across all computational kernels
  int nevt;
#endif

  using mgOnGpu::nw6;

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  inline const fptype& pIparIp4Ievt( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
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
#ifdef __CUDACC__
    const int nevt = blockDim.x * gridDim.x;
#else
    using MG5_sm::nevt;
#endif
    // SOA: allmomenta[npar][np4][ndim]
    return allmomenta[ipar*np4*nevt + ip4*nevt + ievt]; // SOA[ipar][ip4][ievt]
#elif defined MGONGPU_LAYOUT_AOS
    using mgOnGpu::npar;
    // AOS: allmomenta[ndim][npar][np4]
    return allmomenta[ievt*npar*np4 + ipar*np4 + ip4]; // AOS[ievt][ipar][ip4]
#endif
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void imzxxxM0( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
#ifndef __CUDACC__
                 cxtype fi[nw6],
                 const int ievt,
#else
#if defined MGONGPU_WFMEM_LOCAL
                 cxtype fi[nw6],
#else
                 cxtype* fiv,              // output: wavefunctions[5 * 6 * #threads_in_block]
#endif
#endif
                 const int ipar )          // input: particle# out of npar
  {
#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
#if !defined MGONGPU_WFMEM_LOCAL
      const int neib = blockDim.x; // number of events (threads) in block
      const int ieib = threadIdx.x; // index of event (thread) in block
#endif
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "imzxxxM0: ievt=%d ieib=%d\n", ievt, threadIdx.x );
#endif
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
      cxtype& fi0 = fiv[ipar*nw6*neib + 0*neib + ieib];
      cxtype& fi1 = fiv[ipar*nw6*neib + 1*neib + ieib];
      cxtype& fi2 = fiv[ipar*nw6*neib + 2*neib + ieib];
      cxtype& fi3 = fiv[ipar*nw6*neib + 3*neib + ieib];
      cxtype& fi4 = fiv[ipar*nw6*neib + 4*neib + ieib];
      cxtype& fi5 = fiv[ipar*nw6*neib + 5*neib + ieib];
#else
      cxtype& fi0 = fi[0];
      cxtype& fi1 = fi[1];
      cxtype& fi2 = fi[2];
      cxtype& fi3 = fi[3];
      cxtype& fi4 = fi[4];
      cxtype& fi5 = fi[5];
#endif
      fi0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX = PY = 0 and E = -P3 > 0)
      {
        const cxtype chi0 = cxmake( 0, 0 );
        const cxtype chi1 = cxmake( -nhel * sqrt(2 * pvec0), 0 );
        if (nh == 1)
        {
          fi2 = cxmake( 0, 0 );
          fi3 = cxmake( 0, 0 );
          fi4 = chi0;
          fi5 = chi1;
        }
        else
        {
          fi2 = chi1;
          fi3 = chi0;
          fi4 = cxmake( 0, 0 );
          fi5 = cxmake( 0, 0 );
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void ixzxxxM0( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
#ifndef __CUDACC__
                 cxtype fi[nw6],
                 const int ievt,
#else
#if defined MGONGPU_WFMEM_LOCAL
                 cxtype fi[nw6],
#else
                 cxtype* fiv,              // output: wavefunctions[5 * 6 * #threads_in_block]
#endif
#endif
                 const int ipar )          // input: particle# out of npar
  {
#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
#if !defined MGONGPU_WFMEM_LOCAL
      const int neib = blockDim.x; // number of events (threads) in block
      const int ieib = threadIdx.x; // index of event (thread) in block
#endif
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixzxxxM0: ievt=%d ieib=%d\n", ievt, threadIdx.x );
#endif
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
      cxtype& fi0 = fiv[ipar*nw6*neib + 0*neib + ieib];
      cxtype& fi1 = fiv[ipar*nw6*neib + 1*neib + ieib];
      cxtype& fi2 = fiv[ipar*nw6*neib + 2*neib + ieib];
      cxtype& fi3 = fiv[ipar*nw6*neib + 3*neib + ieib];
      cxtype& fi4 = fiv[ipar*nw6*neib + 4*neib + ieib];
      cxtype& fi5 = fiv[ipar*nw6*neib + 5*neib + ieib];
#else
      cxtype& fi0 = fi[0];
      cxtype& fi1 = fi[1];
      cxtype& fi2 = fi[2];
      cxtype& fi3 = fi[3];
      cxtype& fi4 = fi[4];
      cxtype& fi5 = fi[5];
#endif
      fi0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX and PY are not 0)
      {
        const fptype sqp0p3 = sqrt( pvec0 + pvec3 ) * nsf;
        const cxtype chi0 = cxmake( sqp0p3, 0 );
        const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, pvec2 / sqp0p3 );
        if ( nh == 1 )
        {
          fi2 = cxmake( 0, 0 );
          fi3 = cxmake( 0, 0 );
          fi4 = chi0;
          fi5 = chi1;
        }
        else
        {
          fi2 = chi1;
          fi3 = chi0;
          fi4 = cxmake( 0, 0 );
          fi5 = cxmake( 0, 0 );
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void oxzxxxM0( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
#ifndef __CUDACC__
                 cxtype fo[nw6],
                 const int ievt,
#else
#if defined MGONGPU_WFMEM_LOCAL
                 cxtype fo[nw6],
#else
                 cxtype* fov,              // output: wavefunctions[5 * 6 * #threads_in_block]
#endif
#endif
                 const int ipar )          // input: particle# out of npar
  {
#ifndef __CUDACC__
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
#if !defined MGONGPU_WFMEM_LOCAL
      const int neib = blockDim.x; // number of events (threads) in block
      const int ieib = threadIdx.x; // index of event (thread) in block
#endif
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxzxxxM0: ievt=%d ieib=%d\n", ievt, threadIdx.x );
#endif
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
      cxtype& fo0 = fov[ipar*nw6*neib + 0*neib + ieib];
      cxtype& fo1 = fov[ipar*nw6*neib + 1*neib + ieib];
      cxtype& fo2 = fov[ipar*nw6*neib + 2*neib + ieib];
      cxtype& fo3 = fov[ipar*nw6*neib + 3*neib + ieib];
      cxtype& fo4 = fov[ipar*nw6*neib + 4*neib + ieib];
      cxtype& fo5 = fov[ipar*nw6*neib + 5*neib + ieib];
#else
      cxtype& fo0 = fo[0];
      cxtype& fo1 = fo[1];
      cxtype& fo2 = fo[2];
      cxtype& fo3 = fo[3];
      cxtype& fo4 = fo[4];
      cxtype& fo5 = fo[5];
#endif
      fo0 = cxmake( pvec0 * nsf, pvec3 * nsf );
      fo1 = cxmake( pvec1 * nsf, pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // EITHER (Px and Py are not zero)
      // OR (PX = PY = 0 and E = P3 > 0)
      {
        const fptype sqp0p3 = sqrt( pvec0 + pvec3 ) * nsf;
        const cxtype chi0 = cxmake( sqp0p3, 0 );
        const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
        if( nh == 1 )
        {
          fo2 = chi0;
          fo3 = chi1;
          fo4 = cxmake( 0, 0 );
          fo5 = cxmake( 0, 0 );
        }
        else
        {
          fo2 = cxmake( 0, 0 );
          fo3 = cxmake( 0, 0 );
          fo4 = chi1;
          fo5 = chi0;
        }
      }
    }
    // ** END LOOP ON IEVT **
    return;
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void FFV1_0(
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
              const cxtype FAV[],  // input: wavefunctionA[6 * #threads_in_block]
              const cxtype FBV[],  // input: wavefunctionB[6 * #threads_in_block]
              const cxtype VCV[],  // input: wavefunctionC[6 * #threads_in_block]
#else
              const cxtype FA[],   // input wavefunctionA[6]
              const cxtype FB[],   // input wavefunctionB[6]
              const cxtype VC[],   // input wavefunctionC[6]
#endif
              const cxtype COUP,
              cxtype* vertex )    // output
  {
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtype& FA2 = FAV[2*neib + ieib];
    const cxtype& FA3 = FAV[3*neib + ieib];
    const cxtype& FA4 = FAV[4*neib + ieib];
    const cxtype& FA5 = FAV[5*neib + ieib];
    const cxtype& FB2 = FBV[2*neib + ieib];
    const cxtype& FB3 = FBV[3*neib + ieib];
    const cxtype& FB4 = FBV[4*neib + ieib];
    const cxtype& FB5 = FBV[5*neib + ieib];
    const cxtype& VC2 = VCV[2*neib + ieib];
    const cxtype& VC3 = VCV[3*neib + ieib];
    const cxtype& VC4 = VCV[4*neib + ieib];
    const cxtype& VC5 = VCV[5*neib + ieib];
#else
    const cxtype& FA2 = FA[2];
    const cxtype& FA3 = FA[3];
    const cxtype& FA4 = FA[4];
    const cxtype& FA5 = FA[5];
    const cxtype& FB2 = FB[2];
    const cxtype& FB3 = FB[3];
    const cxtype& FB4 = FB[4];
    const cxtype& FB5 = FB[5];
    const cxtype& VC2 = VC[2];
    const cxtype& VC3 = VC[3];
    const cxtype& VC4 = VC[4];
    const cxtype& VC5 = VC[5];
#endif
    const cxtype cI = cxmake( 0, 1 );
    const cxtype TMP4 =
      ( FA2 * ( FB4 * ( VC2 + VC5 ) +
                FB5 * ( VC3 + cI * ( VC4 ) ) ) +
        ( FA3 * ( FB4 * ( VC3 - cI * ( VC4) )
                  + FB5 * ( VC2 - VC5 ) ) +
          ( FA4 * ( FB2 * ( VC2 - VC5 )
                    - FB3 * ( VC3 + cI * ( VC4 ) ) ) +
            FA5 * ( FB2 * ( -VC3 + cI * ( VC4 ) )
                    + FB3 * ( VC2 + VC5 ) )
            )
          )
        );
    ( *vertex ) = COUP * ( -cI ) * TMP4;
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void FFV1P0_3(
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
                const cxtype FAV[],  // input: wavefunctionA[6 * #threads_in_block]
                const cxtype FBV[],  // input: wavefunctionB[6 * #threads_in_block]
#else
                const cxtype FA[],   // input wavefunctionA[6]
                const cxtype FB[],   // input wavefunctionB[6]
#endif
                const cxtype COUP,
                const fptype M3,
                const fptype W3
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
                , cxtype VCV[]         // output: wavefunctionC[6 * #threads_in_block]
#else
                , cxtype VC[]          // output wavefunctionC[6]
#endif
                )
  {
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtype& FA0 = FAV[0*neib + ieib];
    const cxtype& FA1 = FAV[1*neib + ieib];
    const cxtype& FA2 = FAV[2*neib + ieib];
    const cxtype& FA3 = FAV[3*neib + ieib];
    const cxtype& FA4 = FAV[4*neib + ieib];
    const cxtype& FA5 = FAV[5*neib + ieib];
    const cxtype& FB0 = FBV[0*neib + ieib];
    const cxtype& FB1 = FBV[1*neib + ieib];
    const cxtype& FB2 = FBV[2*neib + ieib];
    const cxtype& FB3 = FBV[3*neib + ieib];
    const cxtype& FB4 = FBV[4*neib + ieib];
    const cxtype& FB5 = FBV[5*neib + ieib];
    cxtype& VC0 = VCV[0*neib + ieib];
    cxtype& VC1 = VCV[1*neib + ieib];
    cxtype& VC2 = VCV[2*neib + ieib];
    cxtype& VC3 = VCV[3*neib + ieib];
    cxtype& VC4 = VCV[4*neib + ieib];
    cxtype& VC5 = VCV[5*neib + ieib];
#else
    const cxtype& FA0 = FA[0];
    const cxtype& FA1 = FA[1];
    const cxtype& FA2 = FA[2];
    const cxtype& FA3 = FA[3];
    const cxtype& FA4 = FA[4];
    const cxtype& FA5 = FA[5];
    const cxtype& FB0 = FB[0];
    const cxtype& FB1 = FB[1];
    const cxtype& FB2 = FB[2];
    const cxtype& FB3 = FB[3];
    const cxtype& FB4 = FB[4];
    const cxtype& FB5 = FB[5];
    cxtype& VC0 = VC[0];
    cxtype& VC1 = VC[1];
    cxtype& VC2 = VC[2];
    cxtype& VC3 = VC[3];
    cxtype& VC4 = VC[4];
    cxtype& VC5 = VC[5];
#endif
    const cxtype cI = cxmake( 0, 1 );
    VC0 = + FA0 + FB0;
    VC1 = + FA1 + FB1;
    const fptype PPP0 = -cxreal( VC0 );
    const fptype PPP1 = -cxreal( VC1 );
    const fptype PPP2 = -cximag( VC1 );
    const fptype PPP3 = -cximag( VC0 );
    const cxtype denom =
      COUP / ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) - ( PPP3 * PPP3 ) - M3 * ( M3 - cI * W3 ) );
    VC2 = denom * ( -cI ) * ( FA2 * FB4 + FA3 * FB5 + FA4 * FB2 + FA5 * FB3 );
    VC3 = denom * ( -cI ) * ( -FA2 * FB5 - FA3 * FB4 + FA4 * FB3 + FA5 * FB2 );
    VC4 = denom * ( -cI ) * ( -cI * ( FA2 * FB5 + FA5 * FB2 ) + cI * ( FA3 * FB4 + FA4 * FB3 ) );
    VC5 = denom * ( -cI ) * ( -FA2 * FB4 - FA5 * FB3 + FA3 * FB5 + FA4 * FB2 );
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void FFV2_4_0(
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
                const cxtype FAV[],  // input: wavefunctionA[6 * #threads_in_block]
                const cxtype FBV[],  // input: wavefunctionB[6 * #threads_in_block]
                const cxtype VCV[],  // input: wavefunctionC[6 * #threads_in_block]
#else
                const cxtype FA[],   // input wavefunctionA[6]
                const cxtype FB[],   // input wavefunctionB[6]
                const cxtype VC[],   // input wavefunctionC[6]
#endif
                const cxtype COUP1,
                const cxtype COUP2,
                cxtype* vertex )    // output
  {
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtype& FA2 = FAV[2*neib + ieib];
    const cxtype& FA3 = FAV[3*neib + ieib];
    const cxtype& FA4 = FAV[4*neib + ieib];
    const cxtype& FA5 = FAV[5*neib + ieib];
    const cxtype& FB2 = FBV[2*neib + ieib];
    const cxtype& FB3 = FBV[3*neib + ieib];
    const cxtype& FB4 = FBV[4*neib + ieib];
    const cxtype& FB5 = FBV[5*neib + ieib];
    const cxtype& VC2 = VCV[2*neib + ieib];
    const cxtype& VC3 = VCV[3*neib + ieib];
    const cxtype& VC4 = VCV[4*neib + ieib];
    const cxtype& VC5 = VCV[5*neib + ieib];
#else
    const cxtype& FA2 = FA[2];
    const cxtype& FA3 = FA[3];
    const cxtype& FA4 = FA[4];
    const cxtype& FA5 = FA[5];
    const cxtype& FB2 = FB[2];
    const cxtype& FB3 = FB[3];
    const cxtype& FB4 = FB[4];
    const cxtype& FB5 = FB[5];
    const cxtype& VC2 = VC[2];
    const cxtype& VC3 = VC[3];
    const cxtype& VC4 = VC[4];
    const cxtype& VC5 = VC[5];
#endif
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    const cxtype cI = cxmake( 0, 1 );
    const cxtype TMP2 =
      ( FA4 * ( FB2 * ( VC2 - VC5 ) - FB3 * ( VC3 + cI * ( VC4 ) ) ) +
        FA5 * ( FB2 * ( -VC3 + cI * ( VC4 ) ) + FB3 * ( VC2 + VC5 ) ) );
    const cxtype TMP0 =
      ( FA2 * ( FB4 * ( VC2 + VC5 ) + FB5 * ( VC3 + cI * ( VC4 ) ) ) +
        FA3 * ( FB4 * ( VC3 - cI * ( VC4 ) ) + FB5 * ( VC2 - VC5 ) ) );
    (*vertex) = -fp1 * ( COUP2 * ( +cI * ( TMP0 ) + fp2 * cI * ( TMP2 ) ) + cI * ( TMP0 * COUP1 ) );
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __device__
#endif
  void FFV2_4_3(
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
                const cxtype FAV[],  // input: wavefunctionA[6 * #threads_in_block]
                const cxtype FBV[],  // input: wavefunctionB[6 * #threads_in_block]
#else
                const cxtype FA[],   // input wavefunctionA[6]
                const cxtype FB[],   // input wavefunctionB[6]
#endif
                const cxtype COUP1,
                const cxtype COUP2,
                const fptype M3,
                const fptype W3
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
                , cxtype VCV[]         // output: wavefunctionC[6 * #threads_in_block]
#else
                , cxtype VC[]          // output wavefunctionC[6]
#endif
                )
  {
#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block
    const int ieib = threadIdx.x; // index of event (thread) in block
    const cxtype& FA0 = FAV[0*neib + ieib];
    const cxtype& FA1 = FAV[1*neib + ieib];
    const cxtype& FA2 = FAV[2*neib + ieib];
    const cxtype& FA3 = FAV[3*neib + ieib];
    const cxtype& FA4 = FAV[4*neib + ieib];
    const cxtype& FA5 = FAV[5*neib + ieib];
    const cxtype& FB0 = FBV[0*neib + ieib];
    const cxtype& FB1 = FBV[1*neib + ieib];
    const cxtype& FB2 = FBV[2*neib + ieib];
    const cxtype& FB3 = FBV[3*neib + ieib];
    const cxtype& FB4 = FBV[4*neib + ieib];
    const cxtype& FB5 = FBV[5*neib + ieib];
    cxtype& VC0 = VCV[0*neib + ieib];
    cxtype& VC1 = VCV[1*neib + ieib];
    cxtype& VC2 = VCV[2*neib + ieib];
    cxtype& VC3 = VCV[3*neib + ieib];
    cxtype& VC4 = VCV[4*neib + ieib];
    cxtype& VC5 = VCV[5*neib + ieib];
#else
    const cxtype& FA0 = FA[0];
    const cxtype& FA1 = FA[1];
    const cxtype& FA2 = FA[2];
    const cxtype& FA3 = FA[3];
    const cxtype& FA4 = FA[4];
    const cxtype& FA5 = FA[5];
    const cxtype& FB0 = FB[0];
    const cxtype& FB1 = FB[1];
    const cxtype& FB2 = FB[2];
    const cxtype& FB3 = FB[3];
    const cxtype& FB4 = FB[4];
    const cxtype& FB5 = FB[5];
    cxtype& VC0 = VC[0];
    cxtype& VC1 = VC[1];
    cxtype& VC2 = VC[2];
    cxtype& VC3 = VC[3];
    cxtype& VC4 = VC[4];
    cxtype& VC5 = VC[5];
#endif
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    const cxtype cI = cxmake( 0, 1 );
    fptype OM3 = 0;
    if ( M3 != 0 ) OM3 = fp1 / ( M3 * M3 );
    VC0 = + FA0 + FB0;
    VC1 = + FA1 + FB1;
    const fptype PPP0 = -cxreal( VC0 );
    const fptype PPP1 = -cxreal( VC1 );
    const fptype PPP2 = -cximag( VC1 );
    const fptype PPP3 = -cximag( VC0 );
    const cxtype TMP1 =
      ( FA2 * ( FB4 * ( PPP0 + PPP3 ) + FB5 * ( PPP1 + cI * ( PPP2 ) ) ) +
        FA3 * ( FB4 * ( PPP1 - cI * ( PPP2 ) ) + FB5 * ( PPP0 - PPP3 ) ) );
    const cxtype TMP3 =
      ( FA4 * ( FB2 * ( PPP0 - PPP3 ) - FB3 * ( PPP1 + cI * ( PPP2 ) ) ) +
        FA5 * ( FB2 * ( -PPP1 + cI * ( PPP2 ) ) + FB3 * ( PPP0 + PPP3 ) ) );
    const cxtype denom =
      fp1 / ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) -
              ( PPP3 * PPP3 ) - M3 * (M3 - cI * W3 ) );
    VC2 = denom * ( -fp2 * cI ) *
      ( COUP2 * ( OM3 * ( -fp1/fp2 ) * PPP0 * ( TMP1 + fp2 * ( TMP3 ) )
                  + ( +fp1/fp2 * ( FA2 * FB4 + FA3 * FB5 ) + FA4 * FB2 + FA5 * FB3 ) )
        + fp1/fp2 * ( COUP1 * ( FA2 * FB4 + FA3 * FB5 - PPP0 * OM3 * TMP1 ) ) );
    VC3 = denom * ( -fp2 * cI ) *
      ( COUP2 * ( OM3 * ( -fp1/fp2 ) * PPP1 * ( TMP1 + fp2 * ( TMP3 ) )
                  + ( -fp1/fp2 * ( FA2 * FB5 + FA3 * FB4 ) + FA4 * FB3 + FA5 * FB2 ) )
        - fp1/fp2 * ( COUP1 * ( FA2 * FB5 + FA3 * FB4 + PPP1 * OM3 * TMP1 ) ) );
    VC4 = denom * cI *
      ( COUP2 * ( OM3 * PPP2 * ( TMP1 + fp2 * ( TMP3 ) )
                  + ( +cI * ( FA2 * FB5 ) - cI * ( FA3 * FB4 )
                      - fp2 * cI * ( FA4 * FB3 )
                      + fp2 * cI * ( FA5 * FB2 ) ) )
        + COUP1 * ( +cI * ( FA2 * FB5 ) - cI * ( FA3 * FB4 ) + PPP2 * OM3 * TMP1 ) );
    VC5 = denom * fp2 * cI *
      ( COUP2 * ( OM3 * fp1/fp2 * PPP3 * ( TMP1 + fp2 * ( TMP3 ) ) +
                  ( +fp1/fp2 * ( FA2 * FB4 ) - fp1/fp2 * ( FA3 * FB5 ) - FA4 * FB2 + FA5 * FB3 ) )
        + fp1/fp2 * ( COUP1 * (FA2 * FB4 + PPP3 * OM3 * TMP1 - FA3 * FB5 ) ) );
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

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{
  using mgOnGpu::np4;
  using mgOnGpu::npar;
  const int ncomb = 16; // #helicity combinations is hardcoded for this process (eemumu: ncomb=16)

#ifdef __CUDACC__
  __device__ __constant__ int cHel[ncomb][npar];
  __device__ __constant__ fptype cIPC[6];  // coupling ?
  __device__ __constant__ fptype cIPD[2];
#else
  static int cHel[ncomb][npar];
  static fptype cIPC[6];  // coupling ?
  static fptype cIPD[2];
#endif

#ifdef __CUDACC__
  __device__ unsigned long long sigmakin_itry = 0; // first iteration over nevt events
  __device__ bool sigmakin_goodhel[ncomb] = { false };
#endif

  //--------------------------------------------------------------------------

  using mgOnGpu::nwf;
  using mgOnGpu::nw6;

#ifdef __CUDACC__
#if defined MGONGPU_WFMEM_GLOBAL

  using mgOnGpu::nbpgMAX;
  // Allocate global or shared memory for the wavefunctions of all (external and internal) particles
  // See https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#allocation-persisting-kernel-launches
  __device__ cxtype* dwf[nbpgMAX]; // device dwf[#blocks][5 * 6 * #threads_in_block]

  __global__
  void sigmakin_alloc()
  {
    // Wavefunctions for this block: bwf[5 * 6 * #threads_in_block]
    cxtype*& bwf = dwf[blockIdx.x];

    // Only the first thread in the block does the allocation (we need one allocation per block)
    if ( threadIdx.x == 0 )
    {
      bwf = (cxtype*)malloc( nwf * nw6 * blockDim.x * sizeof(cxtype) ); // cxtype bwf[5 * 6 * #threads_in_block]
      if ( bwf == NULL )
      {
        printf( "ERROR in sigmakin_alloc (block #%4d): malloc failed\n", blockIdx.x );
        assert( bwf != NULL );
      }
      //else printf( "INFO in sigmakin_alloc (block #%4d): malloc successful\n", blockIdx.x );
    }
    __syncthreads();

    // All threads in the block should see the allocation by now
    assert( bwf != NULL );
  }

  __global__
  void sigmakin_free()
  {
    // Only free from one thread!
    // [NB: if this free is missing, cuda-memcheck fails to detect it]
    // [NB: but if free is called twice, cuda-memcheck does detect it]
    cxtype* bwf = dwf[blockIdx.x];
    if ( threadIdx.x == 0 ) free( bwf );
  }

#elif defined MGONGPU_WFMEM_SHARED

  int sigmakin_sharedmem_nbytes( const int ntpb ) // input: #threads per block
  {
    // Wavefunctions for this block: cxtype bwf[5 * 6 * #threads_in_block]
    const int nbytesBwf = nwf * nw6 * ntpb * sizeof(cxtype);
    return nbytesBwf;
  }

#endif
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
#ifdef __CUDACC__
  __device__
#endif
  // ** NB: allmomenta can have three different layouts
  // ASA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
  // SOA: allmomenta[npar][np4][ndim]
  // AOS: allmomenta[ndim][npar][np4]
  void calculate_wavefunctions( int ihel,
                                const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                                fptype &matrix
#ifndef __CUDACC__
                                , const int ievt
#endif
                                )
  {
#ifndef __CUDACC__
    //printf( "calculate_wavefunctions: ievt %d\n", ievt );
#endif
    cxtype amp[2];
#ifdef __CUDACC__
#if !defined MGONGPU_WFMEM_LOCAL
    // eventually move to same AOSOA everywhere, blocks and threads
#if defined MGONGPU_WFMEM_GLOBAL
    const int iblk = blockIdx.x; // index of block in grid
    cxtype* bwf = dwf[iblk];
#elif defined MGONGPU_WFMEM_SHARED
    extern __shared__ cxtype bwf[];
#endif
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][0], -1, bwf, 0 );
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, bwf, 1 );
    MG5_sm::ixzxxxM0( allmomenta, cHel[ihel][2], -1, bwf, 2 );
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][3], +1, bwf, 3 );
#else //local
    cxtype w[nwf][nw6]; // w[5][6]
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][0], -1, w[0], 0 );
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, w[1], 1 );
    MG5_sm::ixzxxxM0( allmomenta, cHel[ihel][2], -1, w[2], 2 );
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][3], +1, w[3], 3 );
#endif
#else // cpp, not cuda
    cxtype w[nwf][nw6]; // w[5][6]
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][0], -1, w[0], ievt, 0 );
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, w[1], ievt, 1 );
    MG5_sm::ixzxxxM0( allmomenta, cHel[ihel][2], -1, w[2], ievt, 2 );
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][3], +1, w[3], ievt, 3 );
#endif

#if defined __CUDACC__ && !defined MGONGPU_WFMEM_LOCAL
    const int neib = blockDim.x; // number of events (threads) in block

    // Diagram 1
    MG5_sm::FFV1P0_3( &(bwf[1*nw6*neib]),
                      &(bwf[0*nw6*neib]),
                      cxmake( cIPC[0], cIPC[1] ),
                      0.,
                      0.,
                      &(bwf[4*nw6*neib]) );

    MG5_sm::FFV1_0( &(bwf[2*nw6*neib]),
                    &(bwf[3*nw6*neib]),
                    &(bwf[4*nw6*neib]),
                    cxmake( cIPC[0], cIPC[1] ),
                    &amp[0] );

    // Diagram 2
    MG5_sm::FFV2_4_3( &(bwf[1*nw6*neib]),
                      &(bwf[0*nw6*neib]),
                      cxmake( cIPC[2], cIPC[3] ),
                      cxmake( cIPC[4], cIPC[5] ),
                      cIPD[0],
                      cIPD[1],
                      &(bwf[4*nw6*neib]) );

    MG5_sm::FFV2_4_0( &(bwf[2*nw6*neib]),
                      &(bwf[3*nw6*neib]),
                      &(bwf[4*nw6*neib]),
                      cxmake( cIPC[2], cIPC[3] ),
                      cxmake( cIPC[4], cIPC[5] ),
                      &amp[1] );
#else

    // Diagram 1
    MG5_sm::FFV1P0_3( w[1], w[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w[4] );
    MG5_sm::FFV1_0( w[2], w[3], w[4], cxmake( cIPC[0], cIPC[1] ), &amp[0] );

    // Diagram 2
    MG5_sm::FFV2_4_3( w[1], w[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w[4] );
    MG5_sm::FFV2_4_0( w[2], w[3], w[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp[1] );

#endif

    const int ncolor = 1;
    cxtype jamp[ncolor];

    // The color matrix;
    static const fptype denom[ncolor] = {1};
    static const fptype cf[ncolor][ncolor] = {{1}};

    // Calculate color flows
    jamp[0] = -amp[0] - amp[1];

    // Sum and square the color flows to get the matrix element
    for( int icol = 0; icol < ncolor; icol++ )
    {
      cxtype ztemp = cxmake( 0, 0 );
      for( int jcol = 0; jcol < ncolor; jcol++ )
        ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
      matrix = matrix + cxreal( ztemp * conj( jamp[icol] ) ) / denom[icol];
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
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(int) ) );
#else
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(int) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2*sizeof(fptype) );
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
    static cxtype tIPC[3] = { cxmake( pars->GC_3 ), cxmake( pars->GC_50 ), cxmake( pars->GC_59 ) };
    static fptype tIPD[2] = { (fptype)pars->mdl_MZ, (fptype)pars->mdl_WZ };

#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 3 * sizeof(cxtype ) ) );
    checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(fptype) ) );
#else
    memcpy( cIPC, tIPC, 3 * sizeof(cxtype ) );
    memcpy( cIPD, tIPD, 2 * sizeof(fptype) );
#endif

  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour.

  // ** NB: allmomenta can have three different layouts
  // ASA: allmomenta[npag][npar][np4][nepp] where ndim=npag*nepp
  // SOA: allmomenta[npar][np4][ndim]
  // AOS: allmomenta[ndim][npar][np4]
#ifdef __CUDACC__
  __global__
#endif
  void sigmaKin( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 fptype* output            // output[nevt]
#ifdef __CUDACC__
                 // NB: nevt == ndim=gpublocks*gputhreads in CUDA
#else
                 , const int nevt          // input: #events
#endif
                 )
  {
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows
    const int maxtry = 10;
#ifndef __CUDACC__
    static unsigned long long sigmakin_itry = 0; // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = { false };
#endif

#ifndef __CUDACC__
    MG5_sm::nevt = nevt;
    // ** START LOOP ON IEVT **
    for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "sigmakin: ievt %d\n", ievt );
#endif

      // Denominators: spins, colors and identical particles
      const int nprocesses = 1;
      const int denominators[nprocesses] = {4};

      // Reset the matrix elements
      fptype matrix_element[nprocesses];
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] = 0.;
      }

      fptype melast = matrix_element[0];
      for (int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry>maxtry && !sigmakin_goodhel[ihel] ) continue;
#ifdef __CUDACC__
        calculate_wavefunctions(ihel, allmomenta, matrix_element[0]); // adds ME for ihel to matrix_element[0]
#else
        calculate_wavefunctions(ihel, allmomenta, matrix_element[0], ievt); // adds ME for ihel to matrix_element[0]
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

#ifndef __CUDACC__
      //if ( sigmakin_itry == maxtry )
      //  for (int ihel = 0; ihel < ncomb; ihel++ )
      //    printf( "sigmakin: ihelgood %2d %d\n", ihel, sigmakin_goodhel[ihel] );
      if ( sigmakin_itry <= maxtry )
        sigmakin_itry++;
#else
      if ( sigmakin_itry <= maxtry )
        atomicAdd(&sigmakin_itry, 1);
#endif

    }
    // ** END LOOP ON IEVT **

  }

  //--------------------------------------------------------------------------

}
