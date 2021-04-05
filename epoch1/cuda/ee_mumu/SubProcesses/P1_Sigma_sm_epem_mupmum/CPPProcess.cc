//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"
#include "mgOnGpuVectors.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <iostream>

mgDebugDeclare();

namespace MG5_sm
{
#ifdef __CUDACC__
  __device__
  inline const fptype& pIparIp4Ievt( const fptype* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                     const int ipar,
                                     const int ip4,
                                     const int ievt )
  {
    // mapping for the various schemes (AOSOA, AOS, SOA...)
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    const int neppM = mgOnGpu::neppM; // AOSOA layout: constant at compile-time
    const int ipagM = ievt/neppM; // #eventpage in this iteration
    const int ieppM = ievt%neppM; // #event in the current eventpage in this iteration
    //printf( "%f\n", momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM] );
    return momenta1d[ipagM*npar*np4*neppM + ipar*np4*neppM + ip4*neppM + ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }
#else
  // Return by value: it seems a tiny bit faster than returning a reference (both for scalar and vector), not clear why
  // NB: this assumes that neppV == neppM!
  inline fptype_sv pIparIp4Ipag( const fptype_sv* momenta1d, // input: momenta as AOSOA[npagM][npar][4][neppM]
                                 const int ipar,
                                 const int ip4,
                                 const int ipagM )
  {
    // mapping for the various schemes (AOSOA, AOS, SOA...)
    using mgOnGpu::np4;
    using mgOnGpu::npar;
    //printf( "%f\n", momenta1d[ipagM*npar*np4 + ipar*np4 + ip4] );
    return momenta1d[ipagM*npar*np4 + ipar*np4 + ip4]; // AOSOA[ipagM][ipar][ip4][ieppM]
  }
#endif

  //--------------------------------------------------------------------------

  /*
  __device__
  void ixxxxx( const fptype* allmomenta,
               const fptype fmass,
               const int nhel,
               const int nsf,
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )  // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
#ifndef __CUDACC__
    using std::max;
    using std::min;
#endif
    // const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
    const fptype& pvec1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt);
    const fptype& pvec2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt);
    const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
    cxtype chi[2];
    fptype sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2];
    int ip, im, nh;
    fptype p[4] = {0, pvec1, pvec2, pvec3};
    p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + fmass * fmass);
    fi[0] = cxtype(-p[0] * nsf, -p[3] * nsf);
    fi[1] = cxtype(-p[1] * nsf, -p[2] * nsf);
    nh = nhel * nsf;
    if (fmass != 0.0)
    {
      pp = min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]));
      if (pp == 0.0)
      {
        sqm[0] = sqrt(std::abs(fmass));
        sqm[1] = (fmass < 0) ? - abs(sqm[0]) : abs(sqm[0]);
        ip = (1 + nh)/2;
        im = (1 - nh)/2;
        fi[2] = ip * sqm[ip];
        fi[3] = im * nsf * sqm[ip];
        fi[4] = ip * nsf * sqm[im];
        fi[5] = im * sqm[im];
      }
      else
      {
        sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5;
        sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5;
        omega[0] = sqrt(p[0] + pp);
        omega[1] = fmass/omega[0];
        ip = (1 + nh)/2;
        im = (1 - nh)/2;
        sfomega[0] = sf[0] * omega[ip];
        sfomega[1] = sf[1] * omega[im];
        pp3 = max(pp + p[3], 0.0);
        chi[0] = cxtype(sqrt(pp3 * 0.5/pp), 0);
        if (pp3 == 0.0)
        {
          chi[1] = cxtype(-nh, 0);
        }
        else
        {
          chi[1] =
            cxtype(nh * p[1], p[2])/sqrt(2.0 * pp * pp3);
        }
        fi[2] = sfomega[0] * chi[im];
        fi[3] = sfomega[0] * chi[ip];
        fi[4] = sfomega[1] * chi[im];
        fi[5] = sfomega[1] * chi[ip];
      }
    }
    else
    {
      if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
      {
        sqp0p3 = 0.0;
      }
      else
      {
        sqp0p3 = sqrt(max(p[0] + p[3], 0.0)) * nsf;
      }
      chi[0] = cxtype(sqp0p3, 0.0);
      if (sqp0p3 == 0.0)
      {
        chi[1] = cxtype(-nhel * sqrt(2.0 * p[0]), 0.0);
      }
      else
      {
        chi[1] = cxtype(nh * p[1], p[2])/sqp0p3;
      }
      if (nh == 1)
      {
        fi[2] = cxtype(0.0, 0.0);
        fi[3] = cxtype(0.0, 0.0);
        fi[4] = chi[0];
        fi[5] = chi[1];
      }
      else
      {
        fi[2] = chi[1];
        fi[3] = chi[0];
        fi[4] = cxtype(0.0, 0.0);
        fi[5] = cxtype(0.0, 0.0);
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void ipzxxx( const fptype * allmomenta,
               const int nhel,
               const int nsf,
               cxtype* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar ) // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
    // const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
    const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
    fi[0] = cxtype (-pvec3 * nsf, -pvec3 * nsf);
    fi[1] = cxtype (0., 0.);
    int nh = nhel * nsf;
    cxtype sqp0p3 = cxtype(sqrt(2. * pvec3) * nsf, 0.);
    fi[2] = fi[1];
    if(nh == 1)
    {
      fi[3] = fi[1];
      fi[4] = sqp0p3;
    }
    else
    {
      fi[3] = sqp0p3;
      fi[4] = fi[1];
    }
    fi[5] = fi[1];
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  __device__
  void imzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const short nhel,
               const short nsf,
               cxtype_sv* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "imzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "imzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      cxtype_sv& fi_0 = fi[0];
      cxtype_sv& fi_1 = fi[1];
      cxtype_sv& fi_2 = fi[2];
      cxtype_sv& fi_3 = fi[3];
      cxtype_sv& fi_4 = fi[4];
      cxtype_sv& fi_5 = fi[5];
      fi_0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi_1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const short nh = nhel * nsf;
      // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
      {
#ifdef __CUDACC__
        const cxtype chi0 = cxmake00();
        const cxtype chi1 = cxmake( -nhel * sqrt(2 * pvec0), 0 );
#else
        const cxtype_sv chi0 = cxmake00();
        const cxtype_sv chi1 = cxmaker0( (short)(-nhel) * sqrt(2 * pvec0) );
#endif
        if (nh == 1)
        {
          fi_2 = cxmake00();
          fi_3 = cxmake00();
          fi_4 = chi0;
          fi_5 = chi1;
        }
        else
        {
          fi_2 = chi1;
          fi_3 = chi0;
          fi_4 = cxmake00();
          fi_5 = cxmake00();
        }
      }
    }
    // +++ END LOOP ON IEVT +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void ixzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const short nhel,
               const short nsf,
               cxtype_sv* fi,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "ixzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "ixzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      cxtype_sv& fi_0 = fi[0];
      cxtype_sv& fi_1 = fi[1];
      cxtype_sv& fi_2 = fi[2];
      cxtype_sv& fi_3 = fi[3];
      cxtype_sv& fi_4 = fi[4];
      cxtype_sv& fi_5 = fi[5];
      fi_0 = cxmake( -pvec0 * nsf, -pvec3 * nsf );
      fi_1 = cxmake( -pvec1 * nsf, -pvec2 * nsf );
      const short nh = nhel * nsf;
      // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
      {
#ifdef __CUDACC__
        const fptype sqp0p3 = sqrt( pvec0 + pvec3 ) * nsf;
        const cxtype chi0 = cxmake( sqp0p3, 0 );
        const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, pvec2 / sqp0p3 );
#else
        const fptype_sv sqp0p3 = sqrt( pvec0 + pvec3 ) * nsf;
        const cxtype_sv chi0 = cxmaker0( sqp0p3 );
        const cxtype_sv chi1 = cxmake( nh * pvec1 / sqp0p3, pvec2 / sqp0p3 );
#endif
        if ( nh == 1 )
        {
          fi_2 = cxmake00();
          fi_3 = cxmake00();
          fi_4 = chi0;
          fi_5 = chi1;
        }
        else
        {
          fi_2 = chi1;
          fi_3 = chi0;
          fi_4 = cxmake00();
          fi_5 = cxmake00();
        }
      }
    }
    // +++ END LOOP ON IEVT +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  /*
  __device__
  void vxxxxx( const fptype* allmomenta,
               const fptype vmass,
               const int nhel,
               const int nsv,
               cxtype* vc,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )  // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    fptype hel, hel0, pt, pt2, pp, pzpt, emp, sqh;
    int nsvahl;
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#else
    using std::min;
#endif
    const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt);
    const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt);
    const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt);
    const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
    // fptype p[4] = {0, pvec[0], pvec[1], pvec[2]};
    // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+vmass*vmass);
    sqh = sqrt(0.5);
    hel = fptype(nhel);
    nsvahl = nsv * std::abs(hel);
    pt2 = (p1 * p1) + (p2 * p2);
    pp = min(p0, sqrt(pt2 + (p3 * p3)));
    pt = min(pp, sqrt(pt2));
    vc[0] = cxtype(p0 * nsv, p3 * nsv);
    vc[1] = cxtype(p1 * nsv, p2 * nsv);
    if (vmass != 0.0)
    {
      hel0 = 1.0 - std::abs(hel);
      if (pp == 0.0)
      {
        vc[2] = cxtype(0.0, 0.0);
        vc[3] = cxtype(-hel * sqh, 0.0);
        vc[4] = cxtype(0.0, nsvahl * sqh);
        vc[5] = cxtype(hel0, 0.0);
      }
      else
      {
        emp = p0/(vmass * pp);
        vc[2] = cxtype(hel0 * pp/vmass, 0.0);
        vc[5] =
          cxtype(hel0 * p3 * emp + hel * pt/pp * sqh, 0.0);
        if (pt != 0.0)
        {
          pzpt = p3/(pp * pt) * sqh * hel;
          vc[3] = cxtype(hel0 * p1 * emp - p1 * pzpt,
                         - nsvahl * p2/pt * sqh);
          vc[4] = cxtype(hel0 * p2 * emp - p2 * pzpt,
                         nsvahl * p1/pt * sqh);
        }
        else
        {
          vc[3] = cxtype(-hel * sqh, 0.0);
          vc[4] = cxtype(0.0, nsvahl * (p3 < 0) ? - abs(sqh)
                         : abs(sqh));
        }
      }
    }
    else
    {
      // pp = p0;
      pt = sqrt((p1 * p1) + (p2 * p2));
      vc[2] = cxtype(0.0, 0.0);
      vc[5] = cxtype(hel * pt/p0 * sqh, 0.0);
      if (pt != 0.0)
      {
        pzpt = p3/(p0 * pt) * sqh * hel;
        vc[3] = cxtype(-p1 * pzpt, -nsv * p2/pt * sqh);
        vc[4] = cxtype(-p2 * pzpt, nsv * p1/pt * sqh);
      }
      else
      {
        vc[3] = cxtype(-hel * sqh, 0.0);
        vc[4] =
          cxtype(0.0, nsv * (p3 < 0) ? - abs(sqh) : abs(sqh));
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void sxxxxx( const fptype* allmomenta,
               const fptype smass,
               const int nhel,
               const int nss,
               cxtype sc[3],
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )
  {
    mgDebug( 0, __FUNCTION__ );
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
    const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt);
    const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt);
    const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt);
    const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
    // fptype p[4] = {0, pvec[0], pvec[1], pvec[2]};
    // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+smass*smass);
    sc[2] = cxtype(1.00, 0.00);
    sc[0] = cxtype(p0 * nss, p3 * nss);
    sc[1] = cxtype(p1 * nss, p2 * nss);
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void oxxxxx( const fptype* allmomenta,
               const fptype fmass,
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )  // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
#ifndef __CUDACC__
    using std::min;
    using std::max;
#endif
    cxtype chi[2];
    fptype sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2];
    int nh, ip, im;
    const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt);
    const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt);
    const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt);
    const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
    // fptype p[4] = {0, pvec[0], pvec[1], pvec[2]};
    // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+fmass*fmass);
    fo[0] = cxtype(p0 * nsf, p3 * nsf);
    fo[1] = cxtype(p1 * nsf, p2 * nsf);
    nh = nhel * nsf;
    if (fmass != 0.000)
    {
      pp = min(p0, sqrt((p1 * p1) + (p2 * p2) + (p3 * p3)));
      if (pp == 0.000)
      {
        sqm[0] = sqrt(std::abs(fmass));
        sqm[1] = (fmass < 0) ? - abs(sqm[0]) : abs(sqm[0]);
        ip = -((1 - nh)/2) * nhel;
        im = (1 + nh)/2 * nhel;
        fo[2] = im * sqm[std::abs(ip)];
        fo[3] = ip * nsf * sqm[std::abs(ip)];
        fo[4] = im * nsf * sqm[std::abs(im)];
        fo[5] = ip * sqm[std::abs(im)];
      }
      else
      {
        sf[0] = fptype(1 + nsf + (1 - nsf) * nh) * 0.5;
        sf[1] = fptype(1 + nsf - (1 - nsf) * nh) * 0.5;
        omega[0] = sqrt(p0 + pp);
        omega[1] = fmass/omega[0];
        ip = (1 + nh)/2;
        im = (1 - nh)/2;
        sfomeg[0] = sf[0] * omega[ip];
        sfomeg[1] = sf[1] * omega[im];
        pp3 = max(pp + p3, 0.00);
        chi[0] = cxtype(sqrt(pp3 * 0.5/pp), 0.00);
        if (pp3 == 0.00)
        {
          chi[1] = cxtype(-nh, 0.00);
        }
        else
        {
          chi[1] =
            cxtype(nh * p1, -p2)/sqrt(2.0 * pp * pp3);
        }
        fo[2] = sfomeg[1] * chi[im];
        fo[3] = sfomeg[1] * chi[ip];
        fo[4] = sfomeg[0] * chi[im];
        fo[5] = sfomeg[0] * chi[ip];
      }
    }
    else
    {
      if ((p1 == 0.00) and (p2 == 0.00) and (p3 < 0.00))
      {
        sqp0p3 = 0.00;
      }
      else
      {
        sqp0p3 = sqrt(max(p0 + p3, 0.00)) * nsf;
      }
      chi[0] = cxtype(sqp0p3, 0.00);
      if (sqp0p3 == 0.000)
      {
        chi[1] = cxtype(-nhel, 0.00) * sqrt(2.0 * p0);
      }
      else
      {
        chi[1] = cxtype(nh * p1, -p2)/sqp0p3;
      }
      if (nh == 1)
      {
        fo[2] = chi[0];
        fo[3] = chi[1];
        fo[4] = cxtype(0.00, 0.00);
        fo[5] = cxtype(0.00, 0.00);
      }
      else
      {
        fo[2] = cxtype(0.00, 0.00);
        fo[3] = cxtype(0.00, 0.00);
        fo[4] = chi[1];
        fo[5] = chi[0];
      }
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void opzxxx( const fptype* allmomenta,
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )  // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == +PZ > 0)
#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
      fo[0] = cxtype (pvec3 * nsf, pvec3 * nsf);
      fo[1] = cxtype (0., 0.);
      int nh = nhel * nsf;
      cxtype CSQP0P3 = cxtype (sqrt(2. * pvec3) * nsf, 0.00);
      fo[3] = fo[1];
      fo[4] = fo[1];
      if (nh == 1)
      {
        fo[2] = CSQP0P3;
        fo[5] = fo[1];
      }
      else
      {
        fo[2] = fo[1];
        fo[5] = CSQP0P3;
      }
    }
    // +++ END LOOP ON IEVT +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void omzxxx( const fptype* allmomenta,
               const int nhel,
               const int nsf,
               cxtype* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )  // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
#ifdef __CUDACC__
    const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
    const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
    fo[0] = cxtype (-pvec3 * nsf, pvec3 * nsf);
    fo[1] = cxtype (0., 0.);
    int nh = nhel * nsf;
    cxtype chi = cxtype (-nhel, 0.00) * sqrt(-2.0 * pvec3);
    if(nh == 1)
    {
      fo[2] = fo[1];
      fo[3] = chi;
      fo[4] = fo[1];
      fo[5] = fo[1];
    }
    else
    {
      fo[2] = fo[1];
      fo[3] = fo[1];
      fo[4] = chi;
      fo[5] = chi;
    }
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  __device__
  void oxzxxx( const fptype_sv* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,        // ASSUME fmass==0
               const short nhel,
               const short nsf,
               cxtype_sv* fo,               // output: wavefunction[(nw6==6)]
#ifndef __CUDACC__
               const int ipagV,
#endif
               const int ipar )             // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and ( either (PT > 0) or (PX == PY == 0 and E == +PZ > 0) )
#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x; // index of event (thread) in grid
      //printf( "oxzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
      const fptype& pvec0 = pIparIp4Ievt( allmomenta, ipar, 0, ievt );
      const fptype& pvec1 = pIparIp4Ievt( allmomenta, ipar, 1, ievt );
      const fptype& pvec2 = pIparIp4Ievt( allmomenta, ipar, 2, ievt );
      const fptype& pvec3 = pIparIp4Ievt( allmomenta, ipar, 3, ievt );
#else
      //printf( "oxzxxx: ipagV=%d\n", ipagV );
      const fptype_sv pvec0 = pIparIp4Ipag( allmomenta, ipar, 0, ipagV );
      const fptype_sv pvec1 = pIparIp4Ipag( allmomenta, ipar, 1, ipagV );
      const fptype_sv pvec2 = pIparIp4Ipag( allmomenta, ipar, 2, ipagV );
      const fptype_sv pvec3 = pIparIp4Ipag( allmomenta, ipar, 3, ipagV );
#endif
      cxtype_sv& fo_0 = fo[0];
      cxtype_sv& fo_1 = fo[1];
      cxtype_sv& fo_2 = fo[2];
      cxtype_sv& fo_3 = fo[3];
      cxtype_sv& fo_4 = fo[4];
      cxtype_sv& fo_5 = fo[5];
      fo_0 = cxmake( pvec0 * nsf, pvec3 * nsf ); // this used to cause the General Protection Fault
      fo_1 = cxmake( pvec1 * nsf, pvec2 * nsf );
      const short nh = nhel * nsf;
      // ASSUMPTIONS: (FMASS == 0) and ( either (PT > 0) or (PX == PY == 0 and E == +PZ > 0) )
      {
#ifdef __CUDACC__
        const fptype sqp0p3 = sqrt( pvec0 + pvec3 ) * nsf;
        const cxtype chi0 = cxmake( sqp0p3, 0 );
        const cxtype chi1 = cxmake( nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
#else
        const fptype_sv sqp0p3 = sqrt( pvec0 + pvec3 ) * nsf;
        const cxtype_sv chi0 = cxmaker0( sqp0p3 );
        const cxtype_sv chi1 = cxmake( nh * pvec1 / sqp0p3, -pvec2 / sqp0p3 );
#endif
        if( nh == 1 )
        {
          fo_2 = chi0;
          fo_3 = chi1;
          fo_4 = cxmake00();
          fo_5 = cxmake00();
        }
        else
        {
          fo_2 = cxmake00();
          fo_3 = cxmake00();
          fo_4 = chi1;
          fo_5 = chi0;
        }
      }
    }
    // +++ END LOOP ON IEVT +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void FFV1_0( const cxtype_sv F1[],   // input: wavefunction1[6]
               const cxtype_sv F2[],   // input: wavefunction2[6]
               const cxtype_sv V3[],   // input: wavefunction3[6]
               const cxtype COUP,
               cxtype_sv* vertex )     // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv& F1_2 = F1[2];
    const cxtype_sv& F1_3 = F1[3];
    const cxtype_sv& F1_4 = F1[4];
    const cxtype_sv& F1_5 = F1[5];
    const cxtype_sv& F2_2 = F2[2];
    const cxtype_sv& F2_3 = F2[3];
    const cxtype_sv& F2_4 = F2[4];
    const cxtype_sv& F2_5 = F2[5];
    const cxtype_sv& V3_2 = V3[2];
    const cxtype_sv& V3_3 = V3[3];
    const cxtype_sv& V3_4 = V3[4];
    const cxtype_sv& V3_5 = V3[5];
    const cxtype cI = cxmake( 0, 1 );
    const cxtype_sv TMP4 =
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

  __device__
  void FFV1P0_3( const cxtype_sv F1[],     // input: wavefunction1[6]
                 const cxtype_sv F2[],     // input: wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )          // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv& F1_0 = F1[0];
    const cxtype_sv& F1_1 = F1[1];
    const cxtype_sv& F1_2 = F1[2];
    const cxtype_sv& F1_3 = F1[3];
    const cxtype_sv& F1_4 = F1[4];
    const cxtype_sv& F1_5 = F1[5];
    const cxtype_sv& F2_0 = F2[0];
    const cxtype_sv& F2_1 = F2[1];
    const cxtype_sv& F2_2 = F2[2];
    const cxtype_sv& F2_3 = F2[3];
    const cxtype_sv& F2_4 = F2[4];
    const cxtype_sv& F2_5 = F2[5];
    cxtype_sv& V3_0 = V3[0];
    cxtype_sv& V3_1 = V3[1];
    cxtype_sv& V3_2 = V3[2];
    cxtype_sv& V3_3 = V3[3];
    cxtype_sv& V3_4 = V3[4];
    cxtype_sv& V3_5 = V3[5];
    const cxtype cI = cxmake( 0, 1 );
    V3_0 = + F1_0 + F2_0;
    V3_1 = + F1_1 + F2_1;
    const fptype_sv PPP0 = -cxreal( V3_0 );
    const fptype_sv PPP1 = -cxreal( V3_1 );
    const fptype_sv PPP2 = -cximag( V3_1 );
    const fptype_sv PPP3 = -cximag( V3_0 );
    const cxtype_sv denom =
      COUP / ( ( PPP0 * PPP0 ) - ( PPP1 * PPP1 ) - ( PPP2 * PPP2 ) - ( PPP3 * PPP3 ) - M3 * ( M3 - cI * W3 ) );
    V3_2 = denom * ( -cI ) * ( F1_2 * F2_4 + F1_3 * F2_5 + F1_4 * F2_2 + F1_5 * F2_3 );
    V3_3 = denom * ( -cI ) * ( -F1_2 * F2_5 - F1_3 * F2_4 + F1_4 * F2_3 + F1_5 * F2_2 );
    V3_4 = denom * ( -cI ) * ( -cI * ( F1_2 * F2_5 + F1_5 * F2_2 ) + cI * ( F1_3 * F2_4 + F1_4 * F2_3 ) );
    V3_5 = denom * ( -cI ) * ( -F1_2 * F2_4 - F1_5 * F2_3 + F1_3 * F2_5 + F1_4 * F2_2 );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  /*
  __device__
  void FFV2_0( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype V3[],   // input: wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )     // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    cxtype cI = cxtype(0., 1.);
    cxtype TMP1;
    TMP1 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
            F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    (*vertex) = COUP * - cI * TMP1;
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void FFV2_3( const cxtype F1[],   // input: wavefunction1[6]
               const cxtype F2[],   // input: wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )        // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    cxtype cI = cxtype(0., 1.);
    fptype OM3;
    fptype P3[4];
    cxtype TMP2;
    cxtype denom;
    OM3 = 0.;
    if (M3 != 0.)
      OM3 = 1./(M3 * M3);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
            F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2);
    V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP2);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2] * OM3 * TMP2);
    V3[5] = denom * (-cI) * (-F1[2] * F2[4] - P3[3] * OM3 * TMP2 + F1[3] * F2[5]);
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void FFV4_0( const cxtype F1[],
               const cxtype F2[],
               const cxtype V3[],
               const cxtype COUP,
               cxtype* vertex )
  {
    mgDebug( 0, __FUNCTION__ );
    cxtype cI = cxtype(0., 1.);
    cxtype TMP3;
    cxtype TMP4;
    TMP3 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
            F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    TMP4 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
            F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])));
    (*vertex) = COUP * (-1.) * (+cI * (TMP3) + 2. * cI * (TMP4));
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  /*
  __device__
  void FFV4_3( const cxtype F1[],
               const cxtype F2[],
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )
  {
    mgDebug( 0, __FUNCTION__ );
    cxtype cI = cxtype(0., 1.);
    fptype OM3;
    fptype P3[4];
    cxtype TMP2;
    cxtype TMP5;
    cxtype denom;
    OM3 = 0.;
    if (M3 != 0.)
      OM3 = 1./(M3 * M3);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    TMP5 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
            F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3])));
    TMP2 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
            F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[0] * (TMP2 + 2. * (TMP5)) +
                                  (+1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] * F2[3]));
    V3[3] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[1] * (TMP2 + 2. * (TMP5)) +
                                  (-1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] * F2[2]));
    V3[4] = denom * 2. * cI * (OM3 * 1./2. * P3[2] * (TMP2 + 2. * (TMP5)) +
                               (+1./2. * cI * (F1[2] * F2[5]) - 1./2. * cI * (F1[3] * F2[4]) - cI *
                                (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
    V3[5] = denom * 2. * cI * (OM3 * 1./2. * P3[3] * (TMP2 + 2. * (TMP5)) +
                               (+1./2. * (F1[2] * F2[4]) - 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] + F1[5] * F2[3]));
    mgDebug( 1, __FUNCTION__ );
    return;
  }
  */

  //--------------------------------------------------------------------------

  __device__
  void FFV2_4_0( const cxtype_sv F1[],   // input: wavefunction1[6]
                 const cxtype_sv F2[],   // input: wavefunction2[6]
                 const cxtype_sv V3[],   // input: wavefunction3[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 cxtype_sv* vertex )     // output: amplitude
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv& F1_2 = F1[2];
    const cxtype_sv& F1_3 = F1[3];
    const cxtype_sv& F1_4 = F1[4];
    const cxtype_sv& F1_5 = F1[5];
    const cxtype_sv& F2_2 = F2[2];
    const cxtype_sv& F2_3 = F2[3];
    const cxtype_sv& F2_4 = F2[4];
    const cxtype_sv& F2_5 = F2[5];
    const cxtype_sv& V3_2 = V3[2];
    const cxtype_sv& V3_3 = V3[3];
    const cxtype_sv& V3_4 = V3[4];
    const cxtype_sv& V3_5 = V3[5];
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    const cxtype cI = cxmake( 0, 1 );
    const cxtype_sv TMP2 =
      ( F1_4 * ( F2_2 * ( V3_2 - V3_5 ) - F2_3 * ( V3_3 + cI * ( V3_4 ) ) ) +
        F1_5 * ( F2_2 * ( -V3_3 + cI * ( V3_4 ) ) + F2_3 * ( V3_2 + V3_5 ) ) );
    const cxtype_sv TMP0 =
      ( F1_2 * ( F2_4 * ( V3_2 + V3_5 ) + F2_5 * ( V3_3 + cI * ( V3_4 ) ) ) +
        F1_3 * ( F2_4 * ( V3_3 - cI * ( V3_4 ) ) + F2_5 * ( V3_2 - V3_5 ) ) );
    (*vertex) = -fp1 * ( COUP2 * ( +cI * ( TMP0 ) + fp2 * cI * ( TMP2 ) ) + cI * ( TMP0 * COUP1 ) );
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void FFV2_4_3( const cxtype_sv F1[],   // input: wavefunction1[6]
                 const cxtype_sv F2[],   // input: wavefunction2[6]
                 const cxtype COUP1,
                 const cxtype COUP2,
                 const fptype M3,
                 const fptype W3,
                 cxtype_sv V3[] )        // output: wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    const cxtype_sv& F1_0 = F1[0];
    const cxtype_sv& F1_1 = F1[1];
    const cxtype_sv& F1_2 = F1[2];
    const cxtype_sv& F1_3 = F1[3];
    const cxtype_sv& F1_4 = F1[4];
    const cxtype_sv& F1_5 = F1[5];
    const cxtype_sv& F2_0 = F2[0];
    const cxtype_sv& F2_1 = F2[1];
    const cxtype_sv& F2_2 = F2[2];
    const cxtype_sv& F2_3 = F2[3];
    const cxtype_sv& F2_4 = F2[4];
    const cxtype_sv& F2_5 = F2[5];
    cxtype_sv& V3_0 = V3[0];
    cxtype_sv& V3_1 = V3[1];
    cxtype_sv& V3_2 = V3[2];
    cxtype_sv& V3_3 = V3[3];
    cxtype_sv& V3_4 = V3[4];
    cxtype_sv& V3_5 = V3[5];
    const fptype fp1 = 1;
    const fptype fp2 = 2;
    const cxtype cI = cxmake( 0, 1 );
    fptype OM3 = 0;
    if ( M3 != 0 ) OM3 = fp1 / ( M3 * M3 );
    V3_0 = + F1_0 + F2_0;
    V3_1 = + F1_1 + F2_1;
    const fptype_sv PPP0 = -cxreal( V3_0 );
    const fptype_sv PPP1 = -cxreal( V3_1 );
    const fptype_sv PPP2 = -cximag( V3_1 );
    const fptype_sv PPP3 = -cximag( V3_0 );
    const cxtype_sv TMP1 =
      ( F1_2 * ( F2_4 * ( PPP0 + PPP3 ) + F2_5 * ( PPP1 + cI * ( PPP2 ) ) ) +
        F1_3 * ( F2_4 * ( PPP1 - cI * ( PPP2 ) ) + F2_5 * ( PPP0 - PPP3 ) ) );
    const cxtype_sv TMP3 =
      ( F1_4 * ( F2_2 * ( PPP0 - PPP3 ) - F2_3 * ( PPP1 + cI * ( PPP2 ) ) ) +
        F1_5 * ( F2_2 * ( -PPP1 + cI * ( PPP2 ) ) + F2_3 * ( PPP0 + PPP3 ) ) );
    const cxtype_sv denom =
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

//#include "../../src/HelAmps_sm.cc"

#include <algorithm>
#include <iostream>

#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

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
  using mgOnGpu::npar; // number of particles in total (initial + final)
  using mgOnGpu::ncomb; // number of helicity combinations

  const int nwf = 5; // #wavefunctions: npar (4 external) + 1 (internal, reused for gamma and Z)
  const int nw6 = 6; // dimension of each wavefunction (see KEK 91-11)

#ifdef __CUDACC__
  __device__ __constant__ int cHel[ncomb][npar];
  //__device__ __constant__ fptype cIPC[6];
  //__device__ __constant__ fptype cIPD[2];
  // FIXME: assume process.nprocesses == 1 for the moment
  //__device__ __constant__ int cNGoodHel[1];
  __device__ __constant__ int cNGoodHel;
  __device__ __constant__ int cGoodHel[ncomb];
#else
  static short cHel[ncomb][npar];
  static fptype cIPC[6];
  static fptype cIPD[2];
  // FIXME: assume process.nprocesses == 1 for the moment
  //static int cNGoodHel[1];
  static int cNGoodHel;
  static int cGoodHel[ncomb];
#endif

  //--------------------------------------------------------------------------

  // Evaluate |M|^2 for each subprocess
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum
  // of |M|^2 over helicities for the given event

  __device__
  void calculate_wavefunctions( int ihel,
                                const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                                fptype* allMEs               // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                                , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                                )
  {
    using namespace MG5_sm;
    mgDebug( 0, __FUNCTION__ );
#ifndef __CUDACC__
    //printf( "calculate_wavefunctions: nevt %d\n", nevt );
#endif

#ifdef __CUDACC__
    const fptype cIPC[6] = { 0, -0.30795376724436879, 0, -0.28804415396362731, 0, 0.082309883272248419 };
    const fptype cIPD[2] = { 91.188000000000002, 2.4414039999999999 };
#endif

    // The number of colors
    const int ncolor = 1;

    // The color matrix
    const fptype denom[ncolor] = {1};
    const fptype cf[ncolor][ncolor] = {{1}};

#ifdef __CUDACC__
    // Local variables for the given event (ievt)
    cxtype w[nwf][nw6]; // w[5][6]
    cxtype amp[1]; // was 2
    cxtype jamp[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams
#else
    // Local variables for the given event page (ipagV)
    cxtype_sv w_v[nwf][nw6]; // w_v[5][6]
    cxtype_sv amp_v[1]; // was 2
    cxtype_sv jamp_v[ncolor] = {}; // sum of the invariant amplitudes for all Feynman diagrams
    //for( int icolor = 0; icolor < ncolor; icolor++ ) jamp_v[icolor] = cxmake00();
#endif

#ifndef __CUDACC__
    const int npagV = nevt / neppV;
    // ** START LOOP ON IPAGV **
    // - default(none): no variables are shared by default
    // - shared: as the name says
    // - private: give each thread its own copy, without initialising
    // - firstprivate: give each thread its own copy, and initialise with value from outside
#pragma omp parallel for default(none) shared(allmomenta,allMEs,cf,cHel,cIPC,cIPD,denom,ihel,npagV) private (amp_v,jamp_v,w_v)
    for ( int ipagV = 0; ipagV < npagV; ++ipagV )
#endif
    {
#ifdef __CUDACC__
      oxzxxx( allmomenta, cHel[ihel][0], -1, w[0], 0 );
#else
      oxzxxx( allmomenta, cHel[ihel][0], -1, w_v[0], ipagV, 0 );
#endif

#ifdef __CUDACC__
      imzxxx( allmomenta, cHel[ihel][1], +1, w[1], 1 );
#else
      imzxxx( allmomenta, cHel[ihel][1], +1, w_v[1], ipagV, 1 );
#endif

#ifdef __CUDACC__
      ixzxxx( allmomenta, cHel[ihel][2], -1, w[2], 2 );
#else
      ixzxxx( allmomenta, cHel[ihel][2], -1, w_v[2], ipagV, 2 );
#endif

#ifdef __CUDACC__
      oxzxxx( allmomenta, cHel[ihel][3], +1, w[3], 3 );
#else
      oxzxxx( allmomenta, cHel[ihel][3], +1, w_v[3], ipagV, 3 );
#endif

#ifndef __CUDACC__
      FFV1P0_3( w_v[1], w_v[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w_v[4] );
      // Amplitude(s) for diagram number 1
      FFV1_0( w_v[2], w_v[3], w_v[4], cxmake( cIPC[0], cIPC[1] ), &amp_v[0] );
      //jamp_v[0] -= amp_v[0];
      jamp_v[0] = jamp_v[0] - amp_v[0];

      FFV2_4_3( w_v[1], w_v[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w_v[4] );
      // Amplitude(s) for diagram number 2
      FFV2_4_0( w_v[2], w_v[3], w_v[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp_v[0] );
      //jamp_v[0] -= amp_v[0];
      jamp_v[0] = jamp_v[0] - amp_v[0];

      // ** START LOOP ON IEPPV **
      for ( int ieppV = 0; ieppV < neppV; ++ieppV )
#endif
      {
#ifdef __CUDACC__
        FFV1P0_3( w[1], w[0], cxmake( cIPC[0], cIPC[1] ), 0., 0., w[4] );
        // Amplitude(s) for diagram number 1
        FFV1_0( w[2], w[3], w[4], cxmake( cIPC[0], cIPC[1] ), &amp[0] );
        jamp[0] -= amp[0];

        FFV2_4_3( w[1], w[0], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), cIPD[0], cIPD[1], w[4] );
        // Amplitude(s) for diagram number 2
        FFV2_4_0( w[2], w[3], w[4], cxmake( cIPC[2], cIPC[3] ), cxmake( cIPC[4], cIPC[5] ), &amp[0] );
        jamp[0] -= amp[0];
#endif

#ifdef __CUDACC__
        const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
        const int ievt = idim;
        //printf( "calculate_wavefunctions: ievt %d\n", ievt );
#else
        const int ievt = ipagV*neppV + ieppV;
        //printf( "calculate_wavefunctions: ievt %d\n", ievt );
#endif

#ifndef __CUDACC__
        // Local variables for the given event (ievt)
#ifdef MGONGPU_CPPSIMD
        cxtype jamp[ncolor];
        jamp[0] = jamp_v[0][ieppV];
#else
        cxtype* jamp = jamp_v;
#endif
#endif

        // Sum and square the color flows to get the matrix element
        // (compute |M|^2 by squaring |M|, taking into account colours)
        for( int icol = 0; icol < ncolor; icol++ )
        {
          cxtype ztemp = cxmake( 0, 0 );
          for( int jcol = 0; jcol < ncolor; jcol++ )
            ztemp += cf[icol][jcol] * jamp[jcol];
          // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum of |M|^2 over helicities for the given event
          // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
          allMEs[ievt] += cxreal( ztemp * conj( jamp[icol] ) ) / denom[icol];
        }

        // Store the leading color flows for choice of color
        // for(i=0;i < ncolor; i++)
        // jamp2[0][i] += cxreal( jamp[i]*conj( jamp[i] ) );

        //printf( "calculate_wavefunction: %6d %2d %f\n", ievt, ihel, allMEs[ievt] );
      }
    }

    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess( int /*numiterations*/,
                          int gpublocks,
                          int gputhreads,
                          bool verbose )
    : //m_numiterations( numiterations ),
    gpu_nblocks( gpublocks ),
    gpu_nthreads( gputhreads ),
    dim( gpu_nblocks * gpu_nthreads ),
    m_verbose( verbose )
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
    const short tHel[ncomb][nexternal] =
      { {-1, -1, -1, -1}, {-1, -1, -1, +1}, {-1, -1, +1, -1}, {-1, -1, +1, +1},
        {-1, +1, -1, -1}, {-1, +1, -1, +1}, {-1, +1, +1, -1}, {-1, +1, +1, +1},
        {+1, -1, -1, -1}, {+1, -1, -1, +1}, {+1, -1, +1, -1}, {+1, -1, +1, +1},
        {+1, +1, -1, -1}, {+1, +1, -1, +1}, {+1, +1, +1, -1}, {+1, +1, +1, +1} };
    memcpy( cHel, tHel, ncomb * nexternal * sizeof(short) );
#endif
    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert( sizeof(cxtype) == 2 * sizeof(fptype) );
#ifndef __CUDACC__
    // SANITY CHECK: momenta AOSOA uses vectors with the same size as fptype_v
    assert( neppV == mgOnGpu::neppM );
#endif
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
    if (m_verbose)
    {
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

#ifdef __CUDACC__
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

  // Retrieve the compiler that was used to build this module
  const std::string CPPProcess::getCompiler()
  {
    std::stringstream out;
#ifdef __CUDACC__
#if defined __CUDACC_VER_MAJOR__ && defined __CUDACC_VER_MINOR__ && defined __CUDACC_VER_BUILD__
    out << "nvcc " << __CUDACC_VER_MAJOR__ << "." << __CUDACC_VER_MINOR__ << "." << __CUDACC_VER_BUILD__;
#else
    out << "nvcc UNKNOWN";
#endif
#elif defined __clang__
#if defined __clang_major__ && defined __clang_minor__ && defined __clang_patchlevel__
    out << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#else
    out << "gcc UNKNOWKN";
#endif

#else
#if defined __GNUC__ && defined __GNUC_MINOR__ && defined __GNUC_PATCHLEVEL__
    out << "gcc (GCC) " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#else
    out << "gcc UNKNOWKN";
#endif
#endif
    return out.str();
  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,              // output: allMEs[nevt], final |M|^2 averaged over all helicities
                            bool* isGoodHel )            // output: isGoodHel[ncomb] - device array
  {
    const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
    const int ievt = idim;
    // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
    fptype allMEsLast = 0;
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to running sum of |M|^2 over helicities for the given event(s)
      calculate_wavefunctions( ihel, allmomenta, allMEs );
      if ( allMEs[ievt] != allMEsLast )
      {
        //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
        isGoodHel[ihel] = true;
      }
      allMEsLast = allMEs[ievt]; // running sum up to helicity ihel for event ievt
    }
  }
#else
  void sigmaKin_getGoodHel( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                            fptype* allMEs,              // output: allMEs[nevt], final |M|^2 averaged over all helicities
                            bool* isGoodHel              // output: isGoodHel[ncomb] - device array
                            , const int nevt )           // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
  {
    const int maxtry0 = ( neppV > 10 ? neppV : 10 ); // 10, but at least neppV (otherwise the npagV loop does not even start)
    fptype allMEsLast[maxtry0] = { 0 };
    const int maxtry = std::min( maxtry0, nevt ); // 10, but at most nevt (avoid invalid memory access if nevt<maxtry0)
    for ( int ievt = 0; ievt < maxtry; ++ievt )
    {
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_getGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      calculate_wavefunctions( ihel, allmomenta, allMEs, maxtry );
      for ( int ievt = 0; ievt < maxtry; ++ievt )
      {
        // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
        if ( allMEs[ievt] != allMEsLast[ievt] )
        {
          //if ( !isGoodHel[ihel] ) std::cout << "sigmaKin_getGoodHel ihel=" << ihel << " TRUE" << std::endl;
          isGoodHel[ihel] = true;
        }
        allMEsLast[ievt] = allMEs[ievt]; // running sum up to helicity ihel
      }
    }
  }
#endif

  //--------------------------------------------------------------------------

  void sigmaKin_setGoodHel( const bool* isGoodHel ) // input: isGoodHel[ncomb] - host array
  {
    // FIXME: assume process.nprocesses == 1 for the moment
    //int nGoodHel[1] = { 0 };
    int nGoodHel = 0;
    int goodHel[ncomb] = { 0 };
    for ( int ihel = 0; ihel < ncomb; ihel++ )
    {
      //std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ? " true" : " false" ) << std::endl;
      if ( isGoodHel[ihel] )
      {
        goodHel[nGoodHel] = ihel;
        nGoodHel++;
      }
    }
#ifdef __CUDACC__
    checkCuda( cudaMemcpyToSymbol( cNGoodHel, &nGoodHel, sizeof(int) ) );
    checkCuda( cudaMemcpyToSymbol( cGoodHel, goodHel, ncomb*sizeof(int) ) );
#else
    cNGoodHel = nGoodHel;
    for ( int ihel = 0; ihel < ncomb; ihel++ ) cGoodHel[ihel] = goodHel[ihel];
#endif
  }

  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour
  // FIXME: assume process.nprocesses == 1 (eventually: allMEs[nevt] -> allMEs[nevt*nprocesses]?)

  __global__
  void sigmaKin( const fptype_sv* allmomenta, // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                 fptype* allMEs               // output: allMEs[nevt], final |M|^2 averaged over all helicities
#ifndef __CUDACC__
                 , const int nevt             // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                 )
  {
    mgDebugInitialise();
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows
    // start sigmakin_lines

    // Denominators: spins, colors and identical particles
    //const int nprocesses = 1;
    //const int denominators[nprocesses] = { 4 };
    const int denominators = 4;

#ifdef __CUDACC__
    // Remember: in CUDA this is a kernel for one event, in c++ this processes n events
    const int idim = blockDim.x * blockIdx.x + threadIdx.x; // event# == threadid (previously was: tid)
    const int ievt = idim;
    //printf( "sigmakin: ievt %d\n", ievt );
#endif

    // PART 0 - INITIALISATION (before calculate_wavefunctions)

#ifndef __CUDACC__
    for ( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
      // Reset the "matrix elements" - running sums of |M|^2 over helicities for the given event
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] = 0; // all zeros
    }

    // PART 1 - HELICITY LOOP: CALCULATE WAVEFUNCTIONS
    // (in both CUDA and C++, using precomputed good helicities)
    for ( int ighel = 0; ighel < cNGoodHel; ighel++ )
    {
      const int ihel = cGoodHel[ighel];
#ifdef __CUDACC__
      calculate_wavefunctions( ihel, allmomenta, allMEs );
#else
      calculate_wavefunctions( ihel, allmomenta, allMEs, nevt );
#endif
    }

    // PART 2 - FINALISATION (after calculate_wavefunctions)
#ifndef __CUDACC__
    for ( int ievt = 0; ievt < nevt; ++ievt )
#endif
    {
      // Get the final |M|^2 as an average over helicities/colors of running sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf1d7c/Handout_4_2016-UZH.pdf]
      // FIXME: assume process.nprocesses == 1 for the moment (eventually: need a loop over processes here?)
      allMEs[ievt] /= denominators;
    }
    mgDebugFinalise();
  }

  //--------------------------------------------------------------------------

}
