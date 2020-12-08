//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================
#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#endif
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

mgDebugDeclare();

namespace MG5_sm
{

  using mgOnGpu::nw6;

  //--------------------------------------------------------------------------

  inline const fptype& pIparIp4Ievt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    //mapping for the various scheme AOS, OSA, ...
    
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    const int neppM = mgOnGpu::neppM; // ASA layout: constant at compile-time
    fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM]) momenta1d; // cast to multiD array pointer (AOSOA)
    const int ipagM = ievt/neppM; // #eventpage in this iteration
    const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
    //return allmomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
    return momenta[ipagM][ipar][ip4][ieppM];
  }

  //--------------------------------------------------------------------------

  void imzxxxM0( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
                 cxtype fis[nw6],
#ifndef SYCL_LANGUAGE_VERSION
                 const int ievt,
#endif
                 const int ipar,
                 sycl::nd_item<3> item_ct1) // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
#ifndef SYCL_LANGUAGE_VERSION
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef SYCL_LANGUAGE_VERSION
      const int ievt =
          item_ct1.get_local_range().get(2) * item_ct1.get_group(2) +
          item_ct1.get_local_id(2); // index of event (thread) in grid
      //printf( "imzxxxM0: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      cxtype& fi_0 = fis[0];
      cxtype& fi_1 = fis[1];
      cxtype& fi_2 = fis[2];
      cxtype& fi_3 = fis[3];
      cxtype& fi_4 = fis[4];
      cxtype& fi_5 = fis[5];
      fi_0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi_1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX = PY = 0 and E = -P3 > 0)
      {
        const cxtype chi0 = cxmake( 0, 0 );
        const cxtype chi1 = cxmake(-nhel * sycl::sqrt(2 * pvec0), 0);
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
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  void ixzxxxM0( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
                 cxtype fis[nw6],          // output: wavefunction[6]
#ifndef SYCL_LANGUAGE_VERSION
                 const int ievt,
#endif
                 const int ipar,
                 sycl::nd_item<3> item_ct1) // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
#ifndef SYCL_LANGUAGE_VERSION
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef SYCL_LANGUAGE_VERSION
      const int ievt =
          item_ct1.get_local_range().get(2) * item_ct1.get_group(2) +
          item_ct1.get_local_id(2); // index of event (thread) in grid
      //printf( "ixzxxxM0: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      cxtype& fi_0 = fis[0];
      cxtype& fi_1 = fis[1];
      cxtype& fi_2 = fis[2];
      cxtype& fi_3 = fis[3];
      cxtype& fi_4 = fis[4];
      cxtype& fi_5 = fis[5];
      fi_0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi_1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // (PX and PY are not 0)
      {
        const fptype sqp0p3 = sycl::sqrt(pvec0 + pvec3) * nsf;
        const cxtype chi0 = cxmake( sqp0p3, 0 );
        const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, pvec2 / sqp0p3 );
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
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------
  SYCL_EXTERNAL
  void oxzxxxM0( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
                 //const fptype fmass,
                 const int nhel,
                 const int nsf,
                 cxtype fos[nw6],          // output: wavefunction[6]
#ifndef SYCL_LANGUAGE_VERSION
                 const int ievt,
#endif
                 const int ipar,
                 sycl::nd_item<3> item_ct1) // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
#ifndef SYCL_LANGUAGE_VERSION
    // ** START LOOP ON IEVT **
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef SYCL_LANGUAGE_VERSION
      const int ievt =
          item_ct1.get_local_range().get(2) * item_ct1.get_group(2) +
          item_ct1.get_local_id(2); // index of event (thread) in grid
      //printf( "oxzxxxM0: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
      cxtype& fo_0 = fos[0];
      cxtype& fo_1 = fos[1];
      cxtype& fo_2 = fos[2];
      cxtype& fo_3 = fos[3];
      cxtype& fo_4 = fos[4];
      cxtype& fo_5 = fos[5];
      fo_0 = cxmake( pvec0 * nsf, pvec3 * nsf );
      fo_1 = cxmake( pvec1 * nsf, pvec2 * nsf );
      const int nh = nhel * nsf;
      // ASSUMPTIONS FMASS = 0 and
      // EITHER (Px and Py are not zero)
      // OR (PX = PY = 0 and E = P3 > 0)
      {
        const fptype sqp0p3 = sycl::sqrt(pvec0 + pvec3) * nsf;
        const cxtype chi0 = cxmake( sqp0p3, 0 );
        const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
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
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  void FFV1_0( const cxtype F1S[],   // input wavefunction1[6]
               const cxtype F2S[],   // input wavefunction2[6]
               const cxtype V3S[],   // input wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )      // output
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype& F1_2 = F1S[2];
    const cxtype& F1_3 = F1S[3];
    const cxtype& F1_4 = F1S[4];
    const cxtype& F1_5 = F1S[5];
    const cxtype& F2_2 = F2S[2];
    const cxtype& F2_3 = F2S[3];
    const cxtype& F2_4 = F2S[4];
    const cxtype& F2_5 = F2S[5];
    const cxtype& V3_2 = V3S[2];
    const cxtype& V3_3 = V3S[3];
    const cxtype& V3_4 = V3S[4];
    const cxtype& V3_5 = V3S[5];
    const cxtype cI = cxmake( 0, 1 );
    const cxtype TMP4 =
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
    ( *vertex ) = COUP * ( -cI ) * TMP4;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  void FFV1P0_3( const cxtype F1S[],   // input wavefunction1[6]
                 const cxtype F2S[],   // input wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3S[] )        // output wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype& F1_0 = F1S[0];
    const cxtype& F1_1 = F1S[1];
    const cxtype& F1_2 = F1S[2];
    const cxtype& F1_3 = F1S[3];
    const cxtype& F1_4 = F1S[4];
    const cxtype& F1_5 = F1S[5];
    const cxtype& F2_0 = F2S[0];
    const cxtype& F2_1 = F2S[1];
    const cxtype& F2_2 = F2S[2];
    const cxtype& F2_3 = F2S[3];
    const cxtype& F2_4 = F2S[4];
    const cxtype& F2_5 = F2S[5];
    cxtype& V3_0 = V3S[0];
    cxtype& V3_1 = V3S[1];
    cxtype& V3_2 = V3S[2];
    cxtype& V3_3 = V3S[3];
    cxtype& V3_4 = V3S[4];
    cxtype& V3_5 = V3S[5];
    const cxtype cI = cxmake( 0, 1 );
    V3_0 = + F1_0 + F2_0;
    V3_1 = + F1_1 + F2_1;
    const fptype PPP0 = -cxreal( V3_0 );
    const fptype PPP1 = -cxreal( V3_1 );
    const fptype PPP2 = -cximag( V3_1 );
    const fptype PPP3 = -cximag( V3_0 );
    const cxtype denom =
      COUP / ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) - ( PPP3 * PPP3 ) - M3 * ( M3 - cI * W3 ) );
    V3_2 = denom * ( -cI ) * ( F1_2 * F2_4 + F1_3 * F2_5 + F1_4 * F2_2 + F1_5 * F2_3 );
    V3_3 = denom * ( -cI ) * ( -F1_2 * F2_5 - F1_3 * F2_4 + F1_4 * F2_3 + F1_5 * F2_2 );
    V3_4 = denom * ( -cI ) * ( -cI * ( F1_2 * F2_5 + F1_5 * F2_2 ) + cI * ( F1_3 * F2_4 + F1_4 * F2_3 ) );
    V3_5 = denom * ( -cI ) * ( -F1_2 * F2_4 - F1_5 * F2_3 + F1_3 * F2_5 + F1_4 * F2_2 );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  void FFV2_4_0( const cxtype F1S[],   // input wavefunction1[6]
                 const cxtype F2S[],   // input wavefunction2[6]
                 const cxtype V3S[],   // input wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype* vertex )      // output
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype& F1_2 = F1S[2];
    const cxtype& F1_3 = F1S[3];
    const cxtype& F1_4 = F1S[4];
    const cxtype& F1_5 = F1S[5];
    const cxtype& F2_2 = F2S[2];
    const cxtype& F2_3 = F2S[3];
    const cxtype& F2_4 = F2S[4];
    const cxtype& F2_5 = F2S[5];
    const cxtype& V3_2 = V3S[2];
    const cxtype& V3_3 = V3S[3];
    const cxtype& V3_4 = V3S[4];
    const cxtype& V3_5 = V3S[5];
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    const cxtype cI = cxmake( 0, 1 );
    const cxtype TMP2 =
      ( F1_4 * ( F2_2 * ( V3_2 - V3_5 ) - F2_3 * ( V3_3 + cI * ( V3_4 ) ) ) +
        F1_5 * ( F2_2 * ( -V3_3 + cI * ( V3_4 ) ) + F2_3 * ( V3_2 + V3_5 ) ) );
    const cxtype TMP0 =
      ( F1_2 * ( F2_4 * ( V3_2 + V3_5 ) + F2_5 * ( V3_3 + cI * ( V3_4 ) ) ) +
        F1_3 * ( F2_4 * ( V3_3 - cI * ( V3_4 ) ) + F2_5 * ( V3_2 - V3_5 ) ) );
    (*vertex) = -fp1 * ( COUP2 * ( +cI * ( TMP0 ) + fp2 * cI * ( TMP2 ) ) + cI * ( TMP0 * COUP1 ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  void FFV2_4_3( const cxtype F1S[],   // input wavefunction1[6]
                 const cxtype F2S[],   // input wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3S[] )        // output wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype& F1_0 = F1S[0];
    const cxtype& F1_1 = F1S[1];
    const cxtype& F1_2 = F1S[2];
    const cxtype& F1_3 = F1S[3];
    const cxtype& F1_4 = F1S[4];
    const cxtype& F1_5 = F1S[5];
    const cxtype& F2_0 = F2S[0];
    const cxtype& F2_1 = F2S[1];
    const cxtype& F2_2 = F2S[2];
    const cxtype& F2_3 = F2S[3];
    const cxtype& F2_4 = F2S[4];
    const cxtype& F2_5 = F2S[5];
    cxtype& V3_0 = V3S[0];
    cxtype& V3_1 = V3S[1];
    cxtype& V3_2 = V3S[2];
    cxtype& V3_3 = V3S[3];
    cxtype& V3_4 = V3S[4];
    cxtype& V3_5 = V3S[5];
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    const cxtype cI = cxmake( 0, 1 );
    fptype OM3 = 0;
    if ( M3 != 0 ) OM3 = fp1 / ( M3 * M3 );
    V3_0 = + F1_0 + F2_0;
    V3_1 = + F1_1 + F2_1;
    const fptype PPP0 = -cxreal( V3_0 );
    const fptype PPP1 = -cxreal( V3_1 );
    const fptype PPP2 = -cximag( V3_1 );
    const fptype PPP3 = -cximag( V3_0 );
    const cxtype TMP1 =
      ( F1_2 * ( F2_4 * ( PPP0 + PPP3 ) + F2_5 * ( PPP1 + cI * ( PPP2 ) ) ) +
        F1_3 * ( F2_4 * ( PPP1 - cI * ( PPP2 ) ) + F2_5 * ( PPP0 - PPP3 ) ) );
    const cxtype TMP3 =
      ( F1_4 * ( F2_2 * ( PPP0 - PPP3 ) - F2_3 * ( PPP1 + cI * ( PPP2 ) ) ) +
        F1_5 * ( F2_2 * ( -PPP1 + cI * ( PPP2 ) ) + F2_3 * ( PPP0 + PPP3 ) ) );
    const cxtype denom =
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
    mgDebug( 1, __FUNCTION__ );
    return;
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
  using mgOnGpu::np4; // 4: the dimension of 4-momenta (E,px,py,pz)
  using mgOnGpu::npar; // 4: #particles in total (external), e+ e- -> mu+ mu-
  using mgOnGpu::ncomb; // 16: #helicity combinations, 2(spin up/down for fermions)**4(npar)

#ifdef __CUDACC__
  __device__ __constant__ int cHel[ncomb][npar];
  //__device__ __constant__ fptype cIPC[6];
  //__device__ __constant__ fptype cIPD[2];
  __device__ __constant__ int cNGoodHel[1];
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static int cHel[ncomb][npar];
  static fptype cIPC[6];
  static fptype cIPD[2];
#endif

  //--------------------------------------------------------------------------

  using mgOnGpu::nwf;
  using mgOnGpu::nw6;

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
  void calculate_wavefunctions( int ihel,
                                const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                                fptype &meHelSum          // input AND output: running sum of |M|^2 over all helicities for this event
                                , sycl::nd_item<3> item_ct1,
                                const sycl::accessor<int, 2, sycl::access::mode::read_write> cHel
#ifndef SYCL_LANGUAGE_VERSION
                                , const int ievt
#endif
                                )
  {
    mgDebug( 0, __FUNCTION__ );
#ifndef SYCL_LANGUAGE_VERSION
    //printf( "calculate_wavefunctions: ievt %d\n", ievt );
#endif

#ifdef SYCL_LANGUAGE_VERSION
    //const int cHel[ncomb][npar] =
    //  { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
    //    {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
    //    {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
    //    {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
    const fptype cIPC[6] = { 0, -0.30795376724436879, 0, -0.28804415396362731, 0, 0.082309883272248419 };
    const fptype cIPD[2] = { 91.188000000000002, 2.4414039999999999 };
#endif

    cxtype amp[2];
    cxtype w[nwf][nw6]; // w[5][6]

#ifdef SYCL_LANGUAGE_VERSION
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][0], -1, w[0], 0, item_ct1);
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, w[1], 1, item_ct1);
    MG5_sm::ixzxxxM0( allmomenta, cHel[ihel][2], -1, w[2], 2, item_ct1);
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][3], +1, w[3], 3, item_ct1);
#else
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][0], -1, w[0], ievt, 0 );
    MG5_sm::imzxxxM0( allmomenta, cHel[ihel][1], +1, w[1], ievt, 1 );
    MG5_sm::ixzxxxM0( allmomenta, cHel[ihel][2], -1, w[2], ievt, 2 );
    MG5_sm::oxzxxxM0( allmomenta, cHel[ihel][3], +1, w[3], ievt, 3 );
#endif

    // Diagram 1
    MG5_sm::FFV1P0_3( w[1], w[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w[4] );
    MG5_sm::FFV1_0( w[2], w[3], w[4], cxmake( cIPC[0], cIPC[1] ), &amp[0] );

    // Diagram 2
    MG5_sm::FFV2_4_3( w[1], w[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w[4] );
    MG5_sm::FFV2_4_0( w[2], w[3], w[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp[1] );

    const int ncolor = 1;
    cxtype jamp[ncolor];

    // The color matrix;
    const fptype denom[ncolor] = {1};
    const fptype cf[ncolor][ncolor] = {{1}};

    // Calculate color flows
    // (compute M as the sum of the invariant amplitudes for all Feynman diagrams)
    jamp[0] = -amp[0] - amp[1];

    // Sum and square the color flows to get the matrix element
    // (compute |M|^2 by squaring |M|, taking into account colours)
    for( int icol = 0; icol < ncolor; icol++ )
    {
      cxtype ztemp = cxmake( 0, 0 );
      for( int jcol = 0; jcol < ncolor; jcol++ )
        ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
      meHelSum = meHelSum + cxreal( ztemp * conj( jamp[icol] ) ) / denom[icol];
    }

    // Store the leading color flows for choice of color
    // for(i=0;i < ncolor; i++)
    // jamp2[0][i] += cxreal( jamp[i]*conj( jamp[i] ) );

    mgDebug( 1, __FUNCTION__ );
    return;
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
#ifdef __CUDACC__
    // Helicities for the process - nodim
    const int tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
    checkCuda( cudaMemcpyToSymbol( cHel, tHel, ncomb * nexternal * sizeof(int) ) );
#else
    // Helicities for the process - nodim
    const int tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
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

#ifdef SYCL_LANGUAGE_VERSION
    //const cxtype tIPC[3] = { cxmake( pars->GC_3 ), cxmake( pars->GC_50 ), cxmake( pars->GC_59 ) };
    //const fptype tIPD[2] = { (fptype)pars->mdl_MZ, (fptype)pars->mdl_WZ };
    //checkCuda( cudaMemcpyToSymbol( cIPC, tIPC, 3 * sizeof(cxtype ) ) );
    //checkCuda( cudaMemcpyToSymbol( cIPD, tIPD, 2 * sizeof(fptype) ) );
#else
    const cxtype tIPC[3] = { cxmake( pars->GC_3 ), cxmake( pars->GC_50 ), cxmake( pars->GC_59 ) };
    const fptype tIPD[2] = { (fptype)pars->mdl_MZ, (fptype)pars->mdl_WZ };
    memcpy( cIPC, tIPC, 3 * sizeof(cxtype) );
    memcpy( cIPD, tIPD, 2 * sizeof(fptype) );
#endif

    //std::cout << std::setprecision(17) << "tIPC[0] = " << tIPC[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[1] = " << tIPC[1] << std::endl;
    //std::cout << std::setprecision(17) << "tIPC[2] = " << tIPC[2] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[0] = " << tIPD[0] << std::endl;
    //std::cout << std::setprecision(17) << "tIPD[1] = " << tIPD[1] << std::endl;
  }

  //--------------------------------------------------------------------------

#ifdef SYCL_LANGUAGE_VERSION
  SYCL_EXTERNAL
  void sigmaKin_getGoodHel( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            bool* isGoodHel,          // output: isGoodHel[ncomb] - device array
                            sycl::nd_item<3> item_ct1,
                            const sycl::accessor<int, 2, sycl::access::mode::read_write> cHel )
  {
    const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
    fptype meHelSum[nprocesses] = { 0 }; // all zeros
    fptype meHelSumLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
      calculate_wavefunctions(ihel, allmomenta, meHelSum[0], item_ct1, cHel);
      if ( meHelSum[0]>meHelSumLast ) isGoodHel[ihel] = true;
      meHelSumLast = meHelSum[0];
    }
  }
#endif

  //--------------------------------------------------------------------------

#ifdef SYCL_LANGUAGE_VERSION
  void sigmaKin_setGoodHel( const bool* isGoodHel, int* cNGoodHel, int* cGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel[1] = { 0 };
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        goodHel[nGoodHel[0]] = ihel;
        nGoodHel[0]++;
      }
    }
    memcpy( cNGoodHel, nGoodHel, sizeof(int) );
    memcpy( cGoodHel, goodHel, ncomb*sizeof(int) );
  }
#endif

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  SYCL_EXTERNAL
  void sigmaKin( const fptype* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype* allMEs            // output: allMEs[nevt], final |M|^2 averaged over all helicities
                 , sycl::nd_item<3> item_ct1,
                 const sycl::accessor<int, 2, sycl::access::mode::read_write> cHel, 
                 int *cNGoodHel,
                 int *cGoodHel 
#ifndef SYCL_LANGUAGE_VERSION
                 , const int nevt          // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
  {
    mgDebugInitialise();
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows
#ifndef SYCL_LANGUAGE_VERSION
    const int maxtry = 10;
    static unsigned long long sigmakin_itry = 0; // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = { false };
#endif

#ifndef SYCL_LANGUAGE_VERSION
    // ** START LOOP ON IEVT **
    for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef SYCL_LANGUAGE_VERSION
      const int idim =
          item_ct1.get_local_range().get(2) * item_ct1.get_group(2) +
          item_ct1.get_local_id(2); // event# == threadid (previously was: tid)
      const int ievt = idim;
      //printf( "sigmakin: ievt %d\n", ievt );
#endif

      // Denominators: spins, colors and identical particles
      const int nprocesses = 1; // FIXME: assume process.nprocesses == 1
      const int denominators[nprocesses] = { 4 };

      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      fptype meHelSum[nprocesses] = { 0 }; // all zeros

#ifdef SYCL_LANGUAGE_VERSION
      // CUDA - using precomputed good helicities
      for ( int ighel = 0; ighel < cNGoodHel[0]; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        calculate_wavefunctions(ihel, allmomenta, meHelSum[0], item_ct1, cHel);
      }
#else
      // C++ - compute good helicities within this loop
      fptype meHelSumLast = 0; // check for good helicities
      for ( int ihel = 0; ihel < ncomb; ihel++ )
      {
        if ( sigmakin_itry>maxtry && !sigmakin_goodhel[ihel] ) continue;
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
        calculate_wavefunctions( ihel, allmomenta, meHelSum[0], ievt );
        if ( sigmakin_itry<=maxtry )
        {
          if ( !sigmakin_goodhel[ihel] && meHelSum[0]>meHelSumLast ) sigmakin_goodhel[ihel] = true;
          meHelSumLast = meHelSum[0];
        }
      }
#endif

      // Get the final |M|^2 as an average over helicities/colors of the running sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
      for (int iproc = 0; iproc < nprocesses; ++iproc)
        meHelSum[iproc] /= denominators[iproc];

      // Set the final average |M|^2 for this event in the output array for all events
      for (int iproc = 0; iproc < nprocesses; ++iproc)
        allMEs[iproc*nprocesses + ievt] = meHelSum[iproc];

#ifndef SYCL_LANGUAGE_VERSION
      if ( sigmakin_itry <= maxtry )
        sigmakin_itry++;
      //if ( sigmakin_itry == maxtry )
      //  for (int ihel = 0; ihel < ncomb; ihel++ )
      //    printf( "sigmakin: ihelgood %2d %d\n", ihel, sigmakin_goodhel[ihel] );
#endif

    }
    // ** END LOOP ON IEVT **
    mgDebugFinalise();

  }

  //--------------------------------------------------------------------------

}
