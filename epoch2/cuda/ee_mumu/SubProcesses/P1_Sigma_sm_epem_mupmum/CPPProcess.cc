//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

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

  __device__
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

  /*
  __device__
  void ixxxxx( const fptype* allmomenta,
               const fptype& fmass,
               const int& nhel,
               const int& nsf,
               cxtype fi[6],
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
               const int& nhel,
               const int& nsf,
               cxtype fi[6],
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
  void imzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int& nhel,
               const int& nsf,
               cxtype fi[6],
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PX == PY == 0 and E == -PZ > 0)
#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
      fi[0] = cxtype (pvec3 * nsf, -pvec3 * nsf);
      fi[1] = cxtype (0., 0.);
      int nh = nhel * nsf;
      cxtype chi = cxtype (-nhel * sqrt(-2.0 * pvec3), 0.0);
      fi[3] = fi[1];
      fi[4] = fi[1];
      if (nh == 1)
      {
        fi[2] = fi[1];
        fi[5] = chi;
      }
      else
      {
        fi[2] = chi;
        fi[5] = fi[1];
      }
    }
    // +++ END LOOP ON IEVT +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void ixzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int& nhel,
               const int& nsf,
               cxtype fi[6],             // output: wavefunction[6]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
  {
    mgDebug( 0, __FUNCTION__ );
    // ASSUMPTIONS: (FMASS == 0) and (PT > 0)
#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    //for (int ievt = 0; ievt < nevt; ++ievt)
#endif
    {
#ifdef __CUDACC__
      const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif
      const fptype& pvec0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt);
      const fptype& pvec1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt);
      const fptype& pvec2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt);
      const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
      // fptype p[4] = {(float), (float) pvec[0], (float) pvec[1], (float) pvec[2]};
      // p[0] = sqrtf(p[3] * p[3] + p[1] * p[1] + p[2] * p[2]);
      fi[0] = cxtype (-pvec0 * nsf, -pvec2 * nsf);
      fi[1] = cxtype (-pvec0 * nsf, -pvec1 * nsf);
      int nh = nhel * nsf;
      float sqp0p3 = sqrtf(pvec0 + pvec3) * nsf;
      cxtype chi0 = cxtype (sqp0p3, 0.0);
      cxtype chi1 = cxtype (nh * pvec1/sqp0p3, pvec2/sqp0p3);
      cxtype CZERO = cxtype(0., 0.);
      if (nh == 1)
      {
        fi[2] = CZERO;
        fi[3] = CZERO;
        fi[4] = chi0;
        fi[5] = chi1;
      }
      else
      {
        fi[2] = chi1;
        fi[3] = chi0;
        fi[4] = CZERO;
        fi[5] = CZERO;
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
               const fptype& vmass,
               const int& nhel,
               const int& nsv,
               cxtype vc[6],
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
               const fptype& smass,
               const int& nhel,
               const int& nss,
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
               const fptype& fmass,
               const int& nhel,
               const int& nsf,
               cxtype fo[6],
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

  __device__
  void opzxxx( const fptype* allmomenta,
               const int& nhel,
               const int& nsf,
               cxtype fo[6],
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

  //--------------------------------------------------------------------------

  /*
  __device__
  void omzxxx( const fptype* allmomenta,
               const int& nhel,
               const int& nsf,
               cxtype fo[6],
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
  void oxzxxx( const fptype* allmomenta, // input[(npar=4)*(np4=4)*nevt]
               //const fptype fmass,     // ASSUME fmass==0
               const int& nhel,
               const int& nsf,
               cxtype fo[6],             // output: wavefunction[6]
#ifndef __CUDACC__
               const int ievt,
#endif
               const int ipar )          // input: particle# out of npar
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
      //printf( "oxzxxx: ievt=%d threadId=%d\n", ievt, threadIdx.x );
#endif
      const fptype& p0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt);
      const fptype& p1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt);
      const fptype& p2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt);
      const fptype& p3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt);
      // float p[4] = {0, (float) pvec[0], (float) pvec[1], (float) pvec[2]};
      // p[0] = sqrtf(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]);
      fo[0] = cxtype (p0 * nsf, p3 * nsf);
      fo[1] = cxtype (p1 * nsf, p2 * nsf);
      int nh = nhel * nsf;
      float sqp0p3 = sqrtf(p0 + p3) * nsf;
      cxtype chi0 = cxtype (sqp0p3, 0.00);
      cxtype chi1 = cxtype (nh * p1/sqp0p3, -p2/sqp0p3);
      cxtype zero = cxtype (0.00, 0.00);
      if(nh == 1)
      {
        fo[2] = chi0;
        fo[3] = chi1;
        fo[4] = zero;
        fo[5] = zero;
      }
      else
      {
        fo[2] = zero;
        fo[3] = zero;
        fo[4] = chi1;
        fo[5] = chi0;
      }
    }
    // +++ END LOOP ON IEVT +++
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void FFV1_0( const cxtype F1[],   // input wavefunction1[6]
               const cxtype F2[],   // input wavefunction2[6]
               const cxtype V3[],   // input wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )     // output
  {
    mgDebug( 0, __FUNCTION__ );
    cxtype cI = cxtype(0., 1.);
    cxtype TMP0;
    TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
            (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
             (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
              F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])))));
    (*vertex) = COUP * - cI * TMP0;
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void FFV1P0_3( const cxtype F1[],   // input wavefunction1[6]
                 const cxtype F2[],   // input wavefunction2[6]
                 const cxtype COUP,
                 const fptype M3,
                 const fptype W3,
                 cxtype V3[] )        // output wavefunction3[6]
  {
    mgDebug( 0, __FUNCTION__ );
    cxtype cI = cxtype(0., 1.);
    fptype P3[4];
    cxtype denom;
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] +
                             F1[5] * F2[3]);
    V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] +
                             F1[5] * F2[2]);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3]
                                                                           * F2[4] + F1[4] * F2[3]));
    V3[5] = denom * (-cI) * (-F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] +
                             F1[4] * F2[2]);
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  /*
  __device__
  void FFV2_0( const cxtype F1[],   // input wavefunction1[6]
               const cxtype F2[],   // input wavefunction2[6]
               const cxtype V3[],   // input wavefunction3[6]
               const cxtype COUP,
               cxtype* vertex )     // output
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
  void FFV2_3( const cxtype F1[],   // input wavefunction1[6]
               const cxtype F2[],   // input wavefunction2[6]
               const cxtype COUP,
               const fptype M3,
               const fptype W3,
               cxtype V3[] )        // output wavefunction3[6]
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
  void FFV2_4_0( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype V3[],
                 const cxtype COUP1,
                 const cxtype COUP2,
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
    (*vertex) = (-1.) * (COUP2 * (+cI * (TMP3) + 2. * cI * (TMP4)) + cI * (TMP3 * COUP1));
    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  __device__
  void FFV2_4_3( const cxtype F1[],
                 const cxtype F2[],
                 const cxtype COUP1,
                 const cxtype COUP2,
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
    denom = 1./((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-2. * cI) *
      (COUP2 * (OM3 * - 1./2. * P3[0] * (TMP2 + 2. * (TMP5)) +
                (+1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] * F2[3])) +
       1./2. * (COUP1 * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP2)));
    V3[3] = denom * (-2. * cI) *
      (COUP2 * (OM3 * - 1./2. * P3[1] * (TMP2 + 2. * (TMP5)) +
                (-1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] * F2[2])) -
       1./2. * (COUP1 * (F1[2] * F2[5] + F1[3] * F2[4] + P3[1] * OM3 * TMP2)));
    V3[4] = denom * cI *
      (COUP2 * (OM3 * P3[2] * (TMP2 + 2. * (TMP5)) +
                (+cI * (F1[2] * F2[5]) - cI * (F1[3] * F2[4]) - 2. * cI * (F1[4] * F2[3]) + 2. * cI * (F1[5] * F2[2]))) +
       COUP1 * (+cI * (F1[2] * F2[5]) - cI * (F1[3] * F2[4]) + P3[2] * OM3 * TMP2));
    V3[5] = denom * 2. * cI *
      (COUP2 * (OM3 * 1./2. * P3[3] * (TMP2 + 2. * (TMP5)) +
                (+1./2. * (F1[2] * F2[4]) - 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] + F1[5] * F2[3])) +
       1./2. * (COUP1 * (F1[2] * F2[4] + P3[3] * OM3 * TMP2 - F1[3] * F2[5])));
    mgDebug( 1, __FUNCTION__ );
    return;
  }

}  // end namespace

//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
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

#ifdef __CUDACC__
  __device__ __constant__ int cHel[ncomb][npar];
  __device__ __constant__ fptype cIPC[6];
  __device__ __constant__ fptype cIPD[2];
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
  // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum
  // of |M|^2 over helicities for the given event

  __device__
  void calculate_wavefunctions( int ihel,
                                const fptype * allmomenta,
                                fptype &meHelSum
#ifndef __CUDACC__
                                , const int ievt
#endif
                                )
  {
    using namespace MG5_sm;
    mgDebug( 0, __FUNCTION__ );
    cxtype amp[1];  // was 2
    const int ncolor = 1;
    cxtype jamp[ncolor];
    // Calculate wavefunctions for all processes
    using namespace MG5_sm;
    cxtype w[nwf][nw6];
    for(int i = 0; i < 1; i++ )
    {
      jamp[i] = cxtype(0., 0.);
    }

#ifdef __CUDACC__
    opzxxx(allmomenta, cHel[ihel][0], -1, w[0], 0);
#else
    opzxxx(allmomenta, cHel[ihel][0], -1, w[0], ievt, 0);
#endif


#ifdef __CUDACC__
    imzxxx(allmomenta, cHel[ihel][1], +1, w[1], 1);
#else
    imzxxx(allmomenta, cHel[ihel][1], +1, w[1], ievt, 1);
#endif


#ifdef __CUDACC__
    ixzxxx(allmomenta, cHel[ihel][2], -1, w[2], 2);
#else
    ixzxxx(allmomenta, cHel[ihel][2], -1, w[2], ievt, 2);
#endif


#ifdef __CUDACC__
    oxzxxx(allmomenta, cHel[ihel][3], +1, w[3], 3);
#else
    oxzxxx(allmomenta, cHel[ihel][3], +1, w[3], ievt, 3);
#endif

    FFV1P0_3(w[1], w[0], cxtype(cIPC[0], cIPC[1]), 0., 0., w[4]);
    // Amplitude(s) for diagram number 1
    FFV1_0(w[2], w[3], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]);
    jamp[0] += -amp[0];
    FFV2_4_3(w[1], w[0], cxtype(cIPC[2], cIPC[3]), cxtype(cIPC[4], cIPC[5]),
             cIPD[0], cIPD[1], w[4]);
    // Amplitude(s) for diagram number 2
    FFV2_4_0(w[2], w[3], w[4], cxtype(cIPC[2], cIPC[3]), cxtype(cIPC[4],
                                                                cIPC[5]), &amp[0]);
    jamp[0] += -amp[0];
    // double CPPProcess::matrix_1_epem_mupmum() {



    // Local variables

    // The color matrix;
    static const double denom[ncolor] = {1};
    static const double cf[ncolor][ncolor] = {{1}};


    // Sum and square the color flows to get the matrix element
    for(int icol = 0; icol < ncolor; icol++ )
    {
      cxtype ztemp = cxmake(0, 0);
      for(int jcol = 0; jcol < ncolor; jcol++ )
        ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
      meHelSum = meHelSum + cxreal(ztemp * conj(jamp[icol]))/denom[icol];
    }

    // Store the leading color flows for choice of color
    // for(i=0;i < ncolor; i++)
    // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

    mgDebug( 1, __FUNCTION__ );
    return;
  }

  //--------------------------------------------------------------------------

  CPPProcess::CPPProcess(int numiterations, int gpublocks, int gputhreads,
                         bool verbose, bool debug)
    : m_numiterations(numiterations), gpu_nblocks(gpublocks),
      gpu_nthreads(gputhreads), m_verbose(verbose),
      dim(gpu_nblocks * gpu_nthreads)
  {


    // Helicities for the process - nodim
    static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1, -1, 1},
                                               {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1}, {-1, 1,
                                                                                                                  1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1, 1, -1},
                                               {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1, 1, 1,
                                                                                                             1}};
#ifdef __CUDACC__
    checkCuda(cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int)));
#else
    memcpy(cHel, tHel, ncomb * nexternal * sizeof(int));
#endif

    // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
    assert(sizeof(cxtype) == 2 * sizeof(fptype));
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

    static cxtype tIPC[3] = {pars->GC_3, pars->GC_50, pars->GC_59};
    static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ};

#ifdef __CUDACC__
    checkCuda(cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(cxtype)));
    checkCuda(cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(fptype)));
#else
    memcpy(cIPC, tIPC, 3 * sizeof(cxtype));
    memcpy(cIPD, tIPD, 2 * sizeof(fptype));
#endif

  }

  //--------------------------------------------------------------------------

#ifdef __CUDACC__
  __global__
  void sigmaKin_getGoodHel(const fptype * allmomenta,  // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
                           bool * isGoodHel)  // output: isGoodHel[ncomb] - device array
  {
    const int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
    fptype meHelSum[nprocesses] = {0};  // all zeros
    fptype meHelSumLast = 0;
    for (int ihel = 0; ihel < ncomb; ihel++ )
    {
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running
      // sum of |M|^2 over helicities for the given event
      calculate_wavefunctions(ihel, allmomenta, meHelSum[0]);
      if (meHelSum[0] != meHelSumLast)
      {
        isGoodHel[ihel] = true;
        meHelSumLast = meHelSum[0];
      }
    }
  }
#endif

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
  void sigmaKin_setGoodHel(const bool * isGoodHel)  // input: isGoodHel[ncomb] - host array
  {
    int nGoodHel[1] = {0};
    int goodHel[ncomb] = {0};
    for (int ihel = 0; ihel < ncomb; ihel++ )
    {
      // std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ?
      // " true" : " false" ) << std::endl;
      if (isGoodHel[ihel])
      {
        goodHel[nGoodHel[0]] = ihel;
        nGoodHel[0]++;
      }
    }
    checkCuda(cudaMemcpyToSymbol(cNGoodHel, nGoodHel, sizeof(int)));
    checkCuda(cudaMemcpyToSymbol(cGoodHel, goodHel, ncomb * sizeof(int)));
  }
#endif



  //--------------------------------------------------------------------------
  // Evaluate |M|^2, part independent of incoming flavour.

  __global__ void sigmaKin(const fptype * allmomenta, fptype * allMEs
#ifndef __CUDACC__
                           , const int nevt  // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
                           )
  {
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();

#ifndef __CUDACC__
    const int maxtry = 10;
    static unsigned long long sigmakin_itry = 0;  // first iteration over nevt events
    static bool sigmakin_goodhel[ncomb] = {false};
#endif

    // Reset color flows


    // start sigmakin_lines


    mgDebugInitialise();
    // Set the parameters which change event by event
    // Need to discuss this with Stefan
    // pars->setDependentParameters();
    // pars->setDependentCouplings();
    // Reset color flows

#ifndef __CUDACC__
    // +++ START LOOP ON IEVT +++
    for (int ievt = 0; ievt < nevt; ++ ievt)
#endif
    {
#ifdef __CUDACC__
      const int idim = blockDim.x * blockIdx.x + threadIdx.x;  // event# == threadid (previously was: tid)
      const int ievt = idim;
      // printf( "sigmakin: ievt %d\n", ievt );
#endif

      // Denominators: spins, colors and identical particles
      const int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
      const int denominators[1] = {4};

      // Reset the "matrix elements" - running sums of |M|^2 over helicities for
      // the given event
      fptype meHelSum[nprocesses] = {0};  // all zeros

#ifdef __CUDACC__
      // CUDA - using precomputed good helicities
      for (int ighel = 0; ighel < cNGoodHel[0]; ighel++ )
      {
        const int ihel = cGoodHel[ighel];
        calculate_wavefunctions(ihel, allmomenta, meHelSum[0]);
      }
#else
      // C++ - compute good helicities within this loop
      fptype meHelSumLast = 0;  // check for good helicities
      for (int ihel = 0; ihel < ncomb; ihel++ )
      {
        if (sigmakin_itry > maxtry && !sigmakin_goodhel[ihel])
          continue;
        // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running
        // sum of |M|^2 over helicities for the given event
        calculate_wavefunctions(ihel, allmomenta, meHelSum[0], ievt);
        if (sigmakin_itry <= maxtry)
        {
          if ( !sigmakin_goodhel[ihel] && meHelSum[0] > meHelSumLast)
            sigmakin_goodhel[ihel] = true;
          meHelSumLast = meHelSum[0];
        }
      }
#endif

      // Get the final |M|^2 as an average over helicities/colors of the running
      // sum of |M|^2 over helicities for the given event
      // [NB 'sum over final spins, average over initial spins', eg see
      // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf
      // 1d7c/Handout_4_2016-UZH.pdf]
      for (int iproc = 0; iproc < nprocesses; ++ iproc)
      {
        meHelSum[iproc] /= denominators[iproc];
      }

      // Set the final average |M|^2 for this event in the output array for all
      // events
      for (int iproc = 0; iproc < nprocesses; ++ iproc)
      {
        allMEs[iproc * nprocesses + ievt] = meHelSum[iproc];
      }

#ifndef __CUDACC__
      if (sigmakin_itry <= maxtry)
        sigmakin_itry++;
      // if ( sigmakin_itry == maxtry )
      // for (int ihel = 0; ihel < ncomb; ihel++ )
      // printf( "sigmakin: ihelgood %2d %d\n", ihel, sigmakin_goodhel[ihel] );
#endif
      // +++ END LOOP ON IEVT +++
      mgDebugFinalise();

    }

  }

  //--------------------------------------------------------------------------

}
