//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

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
inline const fptype& pIparIp4Ievt(const fptype * momenta1d,  // input: momenta as AOSOA[npagM][npar][4][neppM]
const int ipar, 
const int ip4, 
const int ievt)
{
  // mapping for the various scheme AOS, OSA, ...

  using mgOnGpu::np4; 
  using mgOnGpu::npar; 
  const int neppM = mgOnGpu::neppM;  // ASA layout: constant at compile-time
  fptype (*momenta)[npar][np4][neppM] = (fptype (*)[npar][np4][neppM])
  momenta1d;  // cast to multiD array pointer (AOSOA)
  const int ipagM = ievt/neppM;  // #eventpage in this iteration
  const int ieppM = ievt%neppM;  // #event in the current eventpage in this iteration
  // return allmomenta[ipagM*npar*np4*neppM + ipar*neppM*np4 + ip4*neppM +
  // ieppM]; // AOSOA[ipagM][ipar][ip4][ieppM]
  return momenta[ipagM][ipar][ip4][ieppM]; 
}

//--------------------------------------------------------------------------

__device__ void ixxxxx(const fptype * allmomenta, const fptype& fmass, const
int& nhel, const int& nsf, 
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  mgDebug(0, __FUNCTION__); 
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif

  const fptype& pvec0 = pIparIp4Ievt(allmomenta, ipar, 0, ievt); 
  const fptype& pvec1 = pIparIp4Ievt(allmomenta, ipar, 1, ievt); 
  const fptype& pvec2 = pIparIp4Ievt(allmomenta, ipar, 2, ievt); 
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  cxtype chi[2]; 
  fptype sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 

  fptype p[4] = {0, pvec0, pvec1, pvec2}; 
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
  //** END LOOP ON IEVT **
  mgDebug(1, __FUNCTION__); 
  return; 
}


__device__ void ipzxxx(const fptype * allmomenta, const int& nhel, const int&
nsf, 
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTION FMASS == 0
  // PX = PY = 0
  // E = P3 (E>0)
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
  fi[3] = (nh == 1) * fi[1] + (nh == -1) * sqp0p3; 
  fi[4] = (nh == 1) * sqp0p3 + (nh == -1) * fi[1]; 
  fi[5] = fi[1]; 
}

__device__ void imzxxx(const fptype * allmomenta, const int& nhel, const int&
nsf, 
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTION FMASS == 0
  // PX = PY = 0
  // E = -P3 (E>0)
  // printf("p3 %f", pvec[2]);
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 
  fi[0] = cxtype (pvec3 * nsf, -pvec3 * nsf); 
  fi[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 
  cxtype chi = cxtype (-nhel * sqrt(-2.0 * pvec3), 0.0); 

  fi[2] = (nh == 1) * fi[1] + (nh == -1) * chi; 
  fi[3] = fi[1]; 
  fi[4] = fi[1]; 
  fi[5] = (nh == 1) * chi + (nh == -1) * fi[1]; 
}

__device__ void ixzxxx(const fptype * allmomenta, const int& nhel, const int&
nsf, 
cxtype fi[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS: FMASS == 0
  // Px and Py are not zero

  // cxtype chi[2];
  // fptype sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2];
  // int ip, im, nh;
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

  fi[2] = (nh == 1) * CZERO + (nh == -1) * chi1; 
  fi[3] = (nh == 1) * CZERO + (nh == -1) * chi0; 
  fi[4] = (nh == 1) * chi0 + (nh == -1) * CZERO; 
  fi[5] = (nh == 1) * chi1 + (nh == -1) * CZERO; 
  return; 
}

__device__ void vxxxxx(const fptype * allmomenta, const fptype& vmass, const
int& nhel, const int& nsv, 
cxtype vc[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  fptype hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 

#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
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
  return; 
}

__device__ void sxxxxx(const fptype * allmomenta, const fptype& smass, const
int& nhel, const int& nss, 
cxtype sc[3], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)
{
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
  return; 
}

__device__ void oxxxxx(const fptype * allmomenta, const fptype& fmass, const
int& nhel, const int& nsf, 
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
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
  return; 
}

__device__ void opzxxx(const fptype * allmomenta, const int& nhel, const int&
nsf, 
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS FMASS =0
  // PX = PY =0
  // E = PZ
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 

  fo[0] = cxtype (pvec3 * nsf, pvec3 * nsf); 
  fo[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 

  cxtype CSQP0P3 = cxtype (sqrt(2. * pvec3) * nsf, 0.00); 

  fo[2] = (nh == 1) * CSQP0P3 + (nh == -1) * fo[1]; 
  fo[3] = fo[1]; 
  fo[4] = fo[1]; 
  fo[5] = (nh == 1) * fo[1] + (nh == -1) * CSQP0P3; 
}


__device__ void omzxxx(const fptype * allmomenta, const int& nhel, const int&
nsf, 
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS FMASS =0
  // PX = PY =0
  // E = -PZ (E>0)
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
#endif  
  const fptype& pvec3 = pIparIp4Ievt(allmomenta, ipar, 3, ievt); 
  fo[0] = cxtype (-pvec3 * nsf, pvec3 * nsf); 
  fo[1] = cxtype (0., 0.); 
  int nh = nhel * nsf; 
  cxtype chi = cxtype (-nhel, 0.00) * sqrt(-2.0 * pvec3); 

  fo[2] = (nh == 1) * fo[1] + (nh == -1) * fo[1]; 
  fo[3] = (nh == 1) * chi + (nh == -1) * fo[1]; ; 
  fo[4] = (nh == 1) * fo[1] + (nh == -1) * chi; 
  fo[5] = (nh == 1) * fo[1] + (nh == -1) * chi; 

  return; 
}

__device__ void oxzxxx(const fptype * allmomenta, const int& nhel, const int&
nsf, 
cxtype fo[6], 
#ifndef __CUDACC__
const int ievt, 
#endif
const int ipar)  // input: particle# out of npar
{
  // ASSUMPTIONS FMASS =0
  // PT > 0
#ifdef __CUDACC__
  const int ievt = blockDim.x * blockIdx.x + threadIdx.x;  // index of event (thread) in grid
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

  fo[2] = (nh == 1) * chi0 + (nh == -1) * zero; 
  fo[3] = (nh == 1) * chi1 + (nh == -1) * zero; 
  fo[4] = (nh == 1) * zero + (nh == -1) * chi1; 
  fo[5] = (nh == 1) * zero + (nh == -1) * chi0; 

  return; 
}
__device__ void FFV1_0(const cxtype F1[], const cxtype F2[], const cxtype V3[], 
const cxtype COUP, cxtype * vertex)
{
  cxtype cI = cxtype(0., 1.); 
  cxtype TMP0; 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) + 
  (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
  F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))))); 
  (*vertex) = COUP * - cI * TMP0; 
}


__device__ void FFV1_1(const cxtype F2[], const cxtype V3[], const cxtype COUP, 
const fptype M1, const fptype W1, cxtype F1[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P1[4]; 
  cxtype denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  F1[2] = denom * cI * (F2[2] * (P1[0] * (-V3[2] + V3[5]) + (P1[1] * (V3[3] - 
  cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) + V3[4]) + P1[3] * (-V3[2] + 
  V3[5])))) + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * (-1.) * 
  (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * 
  (V3[3] + cI * (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * 
  (V3[3] + cI * (V3[4]))))); 
  F1[3] = denom * (-cI) * (F2[2] * (P1[0] * (-V3[3] + cI * (V3[4])) + (P1[1] * 
  (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] * 
  (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * 
  (-1.) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
   * (V3[2] + V3[5])))) + M1 * (F2[4] * (-V3[3] + cI * (V3[4])) + F2[5] * 
  (-V3[2] + V3[5])))); 
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (-V3[3]
   + cI * (V3[4])) + (P1[2] * (-1.) * (+cI * (V3[3]) + V3[4]) - P1[3] * 
  (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * 
  (-V3[2] + V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] * 
  (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (-V3[2] + V3[5]) + F2[3] * 
  (V3[3] + cI * (V3[4]))))); 
  F1[5] = denom * cI * (F2[4] * (P1[0] * (-V3[3] + cI * (V3[4])) + (P1[1] * 
  (V3[2] + V3[5]) + (P1[2] * (-1.) * (+cI * (V3[2] + V3[5])) + P1[3] * 
  (-V3[3] + cI * (V3[4]))))) + (F2[5] * (P1[0] * (-V3[2] + V3[5]) + (P1[1]
   * (V3[3] + cI * (V3[4])) + (P1[2] * (-cI * (V3[3]) + V3[4]) + P1[3] * 
  (-V3[2] + V3[5])))) + M1 * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * 
  (V3[2] + V3[5])))); 
}


__device__ void FFV1_2(const cxtype F1[], const cxtype V3[], const cxtype COUP, 
const fptype M2, const fptype W2, cxtype F2[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P2[4]; 
  cxtype denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - (P2[2] * P2[2]) - (P2[3] * 
  P2[3]) - M2 * (M2 - cI * W2)); 
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-1.) * 
  (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] * 
  (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * 
  (-V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2]) - cI * (V3[5])) + P2[3] * 
  (-V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * 
  (-V3[3] + cI * (V3[4]))))); 
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (-1.) * (V3[3] + cI * (V3[4])) + 
  (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] * 
  (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (-V3[2] + V3[5]) + (P2[1] * 
  (V3[3] - cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) + V3[4]) + P2[3] * 
  (-V3[2] + V3[5])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * 
  (V3[2] + V3[5])))); 
  F2[4] = denom * (-cI) * (F1[4] * (P2[0] * (-V3[2] + V3[5]) + (P2[1] * (V3[3]
   + cI * (V3[4])) + (P2[2] * (-cI * (V3[3]) + V3[4]) + P2[3] * (-V3[2] + 
  V3[5])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * (-1.) * 
  (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
   * (V3[4]))))) + M2 * (F1[2] * (-1.) * (V3[2] + V3[5]) + F1[3] * (-V3[3] + 
  cI * (V3[4]))))); 
  F2[5] = denom * cI * (F1[4] * (P2[0] * (-1.) * (V3[3] + cI * (V3[4])) + 
  (P2[1] * (V3[2] - V3[5]) + (P2[2] * (+cI * (V3[2]) - cI * (V3[5])) + 
  P2[3] * (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + 
  (P2[1] * (-V3[3] + cI * (V3[4])) + (P2[2] * (-1.) * (+cI * (V3[3]) + 
  V3[4]) - P2[3] * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * 
  (V3[4])) + F1[3] * (V3[2] - V3[5])))); 
}


__device__ void FFV1P0_3(const cxtype F1[], const cxtype F2[], const cxtype
COUP, const fptype M3, const fptype W3, cxtype V3[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P3[4]; 
  cxtype denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * 
  P3[3]) - M3 * (M3 - cI * W3)); 
  V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + 
  F1[5] * F2[3]); 
  V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + 
  F1[5] * F2[2]); 
  V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3]
   * F2[4] + F1[4] * F2[3])); 
  V3[5] = denom * (-cI) * (-F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + 
  F1[4] * F2[2]); 
}


__device__ void VVVV4_0(const cxtype V1[], const cxtype V2[], const cxtype
V3[], const cxtype V4[], const cxtype COUP, cxtype * vertex)
{
  cxtype cI = cxtype(0., 1.); 
  cxtype TMP1; 
  cxtype TMP2; 
  cxtype TMP3; 
  cxtype TMP4; 
  TMP3 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP1 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP2 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  (*vertex) = COUP * (-cI * (TMP1 * TMP2) + cI * (TMP3 * TMP4)); 
}


__device__ void VVVV4P0_1(const cxtype V2[], const cxtype V3[], const cxtype
V4[], const cxtype COUP, const fptype M1, const fptype W1, cxtype V1[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P1[4]; 
  cxtype TMP1; 
  cxtype TMP3; 
  cxtype denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP3 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP1 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (-cI * (V3[2] * TMP1) + cI * (V2[2] * TMP3)); 
  V1[3] = denom * (-cI * (V3[3] * TMP1) + cI * (V2[3] * TMP3)); 
  V1[4] = denom * (-cI * (V3[4] * TMP1) + cI * (V2[4] * TMP3)); 
  V1[5] = denom * (-cI * (V3[5] * TMP1) + cI * (V2[5] * TMP3)); 
}


__device__ void VVVV3_0(const cxtype V1[], const cxtype V2[], const cxtype
V3[], const cxtype V4[], const cxtype COUP, cxtype * vertex)
{
  cxtype cI = cxtype(0., 1.); 
  cxtype TMP3; 
  cxtype TMP4; 
  cxtype TMP5; 
  cxtype TMP6; 
  TMP3 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP5 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  (*vertex) = COUP * (-cI * (TMP5 * TMP6) + cI * (TMP3 * TMP4)); 
}


__device__ void VVVV3P0_1(const cxtype V2[], const cxtype V3[], const cxtype
V4[], const cxtype COUP, const fptype M1, const fptype W1, cxtype V1[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P1[4]; 
  cxtype TMP3; 
  cxtype TMP6; 
  cxtype denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP3 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (-cI * (V4[2] * TMP6) + cI * (V2[2] * TMP3)); 
  V1[3] = denom * (-cI * (V4[3] * TMP6) + cI * (V2[3] * TMP3)); 
  V1[4] = denom * (-cI * (V4[4] * TMP6) + cI * (V2[4] * TMP3)); 
  V1[5] = denom * (-cI * (V4[5] * TMP6) + cI * (V2[5] * TMP3)); 
}


__device__ void VVV1_0(const cxtype V1[], const cxtype V2[], const cxtype V3[], 
const cxtype COUP, cxtype * vertex)
{
  cxtype cI = cxtype(0., 1.); 
  fptype P1[4]; 
  fptype P2[4]; 
  fptype P3[4]; 
  cxtype TMP10; 
  cxtype TMP11; 
  cxtype TMP12; 
  cxtype TMP2; 
  cxtype TMP4; 
  cxtype TMP6; 
  cxtype TMP7; 
  cxtype TMP8; 
  cxtype TMP9; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP7 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP8 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP11 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP12 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP4 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP2 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP10 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  (*vertex) = COUP * (TMP2 * (+cI * (TMP9) - cI * (TMP10)) + (TMP4 * (-cI * 
  (TMP7) + cI * (TMP8)) + TMP6 * (-cI * (TMP11) + cI * (TMP12)))); 
}


__device__ void VVV1P0_1(const cxtype V2[], const cxtype V3[], const cxtype
COUP, const fptype M1, const fptype W1, cxtype V1[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P1[4]; 
  fptype P2[4]; 
  fptype P3[4]; 
  cxtype TMP10; 
  cxtype TMP6; 
  cxtype TMP7; 
  cxtype TMP8; 
  cxtype TMP9; 
  cxtype denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP7 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP8 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP9 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP10 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (TMP6 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[2] * (+cI * (TMP9) - cI * (TMP10)))); 
  V1[3] = denom * (TMP6 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[3] * (+cI * (TMP9) - cI * (TMP10)))); 
  V1[4] = denom * (TMP6 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[4] * (+cI * (TMP9) - cI * (TMP10)))); 
  V1[5] = denom * (TMP6 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[5] * (+cI * (TMP9) - cI * (TMP10)))); 
}


__device__ void VVVV1_0(const cxtype V1[], const cxtype V2[], const cxtype
V3[], const cxtype V4[], const cxtype COUP, cxtype * vertex)
{
  cxtype cI = cxtype(0., 1.); 
  cxtype TMP1; 
  cxtype TMP2; 
  cxtype TMP5; 
  cxtype TMP6; 
  TMP1 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP5 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP2 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  (*vertex) = COUP * (-cI * (TMP5 * TMP6) + cI * (TMP1 * TMP2)); 
}


__device__ void VVVV1P0_1(const cxtype V2[], const cxtype V3[], const cxtype
V4[], const cxtype COUP, const fptype M1, const fptype W1, cxtype V1[])
{
  cxtype cI = cxtype(0., 1.); 
  fptype P1[4]; 
  cxtype TMP1; 
  cxtype TMP6; 
  cxtype denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP1 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (-cI * (V4[2] * TMP6) + cI * (V3[2] * TMP1)); 
  V1[3] = denom * (-cI * (V4[3] * TMP6) + cI * (V3[3] * TMP1)); 
  V1[4] = denom * (-cI * (V4[4] * TMP6) + cI * (V3[4] * TMP1)); 
  V1[5] = denom * (-cI * (V4[5] * TMP6) + cI * (V3[5] * TMP1)); 
}



}  // end namespace $(namespace)s_sm





#include <algorithm>
#include <iostream>
#include "mgOnGpuTypes.h"
#include "mgOnGpuConfig.h"

#include "CPPProcess.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ g g WEIGHTED<=4 @1

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{

using mgOnGpu::np4;  // 4: the dimension of 4-momenta (E,px,py,pz)
using mgOnGpu::npar;  // number of particles in total (initial + final)
using mgOnGpu::ncomb;  // number of helicity combinations

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

__device__ void calculate_wavefunctions(int ihel, const fptype * allmomenta,
    fptype &meHelSum
#ifndef __CUDACC__
, const int ievt
#endif
)
{
  using namespace MG5_sm; 
  mgDebug(0, __FUNCTION__); 
  cxtype amp[1];  // was 159
  const int ncolor = 24; 
  cxtype jamp[ncolor]; 
  // Calculate wavefunctions for all processes
  using namespace MG5_sm; 
  cxtype w[nwf][nw6]; 
  for(int i = 0; i < 24; i++ )
  {
    jamp[i] = cxtype(0., 0.); 
  }

#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][0], -1, w[0], 0); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][0], -1, w[0], ievt, 0); 
#endif 


#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][1], -1, w[1], 1); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][1], -1, w[1], ievt, 1); 
#endif 


#ifdef __CUDACC__
  oxxxxx(allmomenta, cIPD[0], cHel[ihel][2], +1, w[2], 2); 
#else
  oxxxxx(allmomenta, cIPD[0], cHel[ihel][2], +1, w[2], ievt, 2); 
#endif 


#ifdef __CUDACC__
  ixxxxx(allmomenta, cIPD[0], cHel[ihel][3], -1, w[3], 3); 
#else
  ixxxxx(allmomenta, cIPD[0], cHel[ihel][3], -1, w[3], ievt, 3); 
#endif 


#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][4], +1, w[4], 4); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][4], +1, w[4], ievt, 4); 
#endif 


#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][5], +1, w[5], 5); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][5], +1, w[5], ievt, 5); 
#endif 

  VVV1P0_1(w[0], w[1], cxtype(cIPC[0], cIPC[1]), 0., 0., w[6]); 
  FFV1P0_3(w[3], w[2], cxtype(cIPC[2], cIPC[3]), 0., 0., w[7]); 
  // Amplitude(s) for diagram number 1
  VVVV1_0(w[6], w[7], w[4], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[6], w[7], w[4], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[6], w[7], w[4], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[6], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[8]); 
  // Amplitude(s) for diagram number 2
  VVV1_0(w[7], w[5], w[8], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[6], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[9]); 
  // Amplitude(s) for diagram number 3
  VVV1_0(w[7], w[4], w[9], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[4], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[10]); 
  // Amplitude(s) for diagram number 4
  VVV1_0(w[6], w[7], w[10], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[11]); 
  FFV1_2(w[3], w[6], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  // Amplitude(s) for diagram number 5
  FFV1_0(w[12], w[11], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 6
  FFV1_0(w[3], w[11], w[9], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +amp[0]; 
  jamp[14] += -amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[17] += +amp[0]; 
  FFV1_2(w[3], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[13]); 
  // Amplitude(s) for diagram number 7
  FFV1_0(w[13], w[11], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[14]); 
  // Amplitude(s) for diagram number 8
  FFV1_0(w[12], w[14], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 9
  FFV1_0(w[3], w[14], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +amp[0]; 
  jamp[20] += -amp[0]; 
  jamp[22] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_2(w[3], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[15]); 
  // Amplitude(s) for diagram number 10
  FFV1_0(w[15], w[14], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[6], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[16]); 
  // Amplitude(s) for diagram number 11
  FFV1_0(w[15], w[16], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 12
  FFV1_0(w[15], w[2], w[9], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[18] += -amp[0]; 
  jamp[20] += +amp[0]; 
  // Amplitude(s) for diagram number 13
  FFV1_0(w[13], w[16], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 14
  FFV1_0(w[13], w[2], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[6] += -amp[0]; 
  jamp[12] += -amp[0]; 
  jamp[14] += +amp[0]; 
  // Amplitude(s) for diagram number 15
  FFV1_0(w[3], w[16], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[6] += -amp[0]; 
  jamp[7] += +amp[0]; 
  // Amplitude(s) for diagram number 16
  FFV1_0(w[12], w[2], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += +amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[22] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_1(w[2], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  FFV1_2(w[3], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[16]); 
  FFV1_1(w[12], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[8]); 
  // Amplitude(s) for diagram number 17
  FFV1_0(w[16], w[8], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += -amp[0]; 
  FFV1_1(w[12], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[9]); 
  // Amplitude(s) for diagram number 18
  FFV1_0(w[16], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += -amp[0]; 
  // Amplitude(s) for diagram number 19
  FFV1_0(w[16], w[12], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[1], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[6]); 
  FFV1P0_3(w[3], w[12], cxtype(cIPC[2], cIPC[3]), 0., 0., w[17]); 
  // Amplitude(s) for diagram number 20
  VVV1_0(w[6], w[5], w[17], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[4] += -amp[0]; 
  jamp[5] += +amp[0]; 
  // Amplitude(s) for diagram number 21
  FFV1_0(w[3], w[9], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 22
  FFV1_0(w[13], w[12], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[1], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[18]); 
  // Amplitude(s) for diagram number 23
  VVV1_0(w[18], w[4], w[17], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[3] += +amp[0]; 
  jamp[4] += -amp[0]; 
  // Amplitude(s) for diagram number 24
  FFV1_0(w[3], w[8], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 25
  FFV1_0(w[15], w[12], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[12], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[19]); 
  // Amplitude(s) for diagram number 26
  FFV1_0(w[15], w[19], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += -amp[0]; 
  // Amplitude(s) for diagram number 27
  FFV1_0(w[15], w[9], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += -amp[0]; 
  // Amplitude(s) for diagram number 28
  FFV1_0(w[13], w[19], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += -amp[0]; 
  // Amplitude(s) for diagram number 29
  FFV1_0(w[13], w[8], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += -amp[0]; 
  // Amplitude(s) for diagram number 30
  FFV1_0(w[3], w[19], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 31
  VVV1_0(w[1], w[10], w[17], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[5] += +amp[0]; 
  VVVV1P0_1(w[1], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[17]); 
  VVVV3P0_1(w[1], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[19]); 
  VVVV4P0_1(w[1], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[8]); 
  // Amplitude(s) for diagram number 32
  FFV1_0(w[3], w[12], w[17], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[5] += +amp[0]; 
  FFV1_0(w[3], w[12], w[19], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[4] += +amp[0]; 
  FFV1_0(w[3], w[12], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[4] += +amp[0]; 
  jamp[5] += -amp[0]; 
  FFV1_2(w[3], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  FFV1_1(w[2], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[9]); 
  FFV1_2(w[12], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[20]); 
  // Amplitude(s) for diagram number 33
  FFV1_0(w[20], w[9], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[11] += -amp[0]; 
  FFV1_2(w[12], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[21]); 
  // Amplitude(s) for diagram number 34
  FFV1_0(w[21], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += -amp[0]; 
  // Amplitude(s) for diagram number 35
  FFV1_0(w[12], w[9], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  FFV1P0_3(w[12], w[2], cxtype(cIPC[2], cIPC[3]), 0., 0., w[22]); 
  // Amplitude(s) for diagram number 36
  VVV1_0(w[6], w[5], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[9] += +amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[23] += +amp[0]; 
  // Amplitude(s) for diagram number 37
  FFV1_0(w[21], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 38
  FFV1_0(w[12], w[14], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 39
  VVV1_0(w[18], w[4], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[11] += +amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[17] += +amp[0]; 
  jamp[21] += -amp[0]; 
  // Amplitude(s) for diagram number 40
  FFV1_0(w[20], w[2], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 41
  FFV1_0(w[12], w[11], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  FFV1_2(w[12], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[23]); 
  // Amplitude(s) for diagram number 42
  FFV1_0(w[23], w[11], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[17] += -amp[0]; 
  // Amplitude(s) for diagram number 43
  FFV1_0(w[21], w[11], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[15] += -amp[0]; 
  // Amplitude(s) for diagram number 44
  FFV1_0(w[23], w[14], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[23] += -amp[0]; 
  // Amplitude(s) for diagram number 45
  FFV1_0(w[20], w[14], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[21] += -amp[0]; 
  // Amplitude(s) for diagram number 46
  FFV1_0(w[23], w[2], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 47
  VVV1_0(w[1], w[10], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[9] += +amp[0]; 
  jamp[11] += -amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[23] += +amp[0]; 
  // Amplitude(s) for diagram number 48
  FFV1_0(w[12], w[2], w[17], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += +amp[0]; 
  jamp[11] += -amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_0(w[12], w[2], w[19], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[11] += -amp[0]; 
  jamp[15] += +amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[21] += +amp[0]; 
  FFV1_0(w[12], w[2], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += -amp[0]; 
  jamp[15] += +amp[0]; 
  jamp[21] += +amp[0]; 
  jamp[23] += -amp[0]; 
  VVV1P0_1(w[0], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[12]); 
  FFV1_2(w[3], w[12], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[22]); 
  // Amplitude(s) for diagram number 49
  FFV1_0(w[22], w[9], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[12], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[23]); 
  // Amplitude(s) for diagram number 50
  FFV1_0(w[3], w[9], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[11] += +amp[0]; 
  // Amplitude(s) for diagram number 51
  FFV1_0(w[13], w[9], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[12], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[20]); 
  // Amplitude(s) for diagram number 52
  FFV1_0(w[16], w[20], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 53
  FFV1_0(w[16], w[2], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[19] += -amp[0]; 
  jamp[22] += +amp[0]; 
  // Amplitude(s) for diagram number 54
  FFV1_0(w[16], w[14], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 55
  FFV1_0(w[3], w[20], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[12] += -amp[0]; 
  jamp[13] += +amp[0]; 
  // Amplitude(s) for diagram number 56
  FFV1_0(w[22], w[2], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += +amp[0]; 
  jamp[11] += -amp[0]; 
  jamp[20] += -amp[0]; 
  jamp[21] += +amp[0]; 
  // Amplitude(s) for diagram number 57
  VVV1_0(w[12], w[18], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 58
  VVVV1_0(w[12], w[1], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[12], w[1], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[12], w[1], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[12], w[1], cxtype(cIPC[0], cIPC[1]), 0., 0., w[21]); 
  // Amplitude(s) for diagram number 59
  VVV1_0(w[7], w[5], w[21], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 60
  VVV1_0(w[1], w[7], w[23], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 61
  FFV1_0(w[3], w[14], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += +amp[0]; 
  jamp[20] += -amp[0]; 
  jamp[21] += +amp[0]; 
  jamp[22] += -amp[0]; 
  // Amplitude(s) for diagram number 62
  FFV1_0(w[22], w[14], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 63
  FFV1_0(w[13], w[2], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +amp[0]; 
  jamp[6] += -amp[0]; 
  jamp[8] += +amp[0]; 
  jamp[12] += -amp[0]; 
  // Amplitude(s) for diagram number 64
  FFV1_0(w[13], w[20], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[20]); 
  FFV1_2(w[3], w[20], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[21]); 
  // Amplitude(s) for diagram number 65
  FFV1_0(w[21], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[20], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[22]); 
  // Amplitude(s) for diagram number 66
  FFV1_0(w[3], w[9], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += +amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[9] += +amp[0]; 
  jamp[10] += -amp[0]; 
  // Amplitude(s) for diagram number 67
  FFV1_0(w[15], w[9], w[20], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[20], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[23]); 
  // Amplitude(s) for diagram number 68
  FFV1_0(w[16], w[23], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 69
  FFV1_0(w[16], w[2], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[16] += +amp[0]; 
  jamp[19] += -amp[0]; 
  // Amplitude(s) for diagram number 70
  FFV1_0(w[16], w[11], w[20], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 71
  FFV1_0(w[3], w[23], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +amp[0]; 
  jamp[5] += -amp[0]; 
  jamp[18] += -amp[0]; 
  jamp[19] += +amp[0]; 
  // Amplitude(s) for diagram number 72
  FFV1_0(w[21], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += +amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[14] += -amp[0]; 
  jamp[15] += +amp[0]; 
  // Amplitude(s) for diagram number 73
  VVV1_0(w[20], w[6], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 74
  VVVV1_0(w[20], w[1], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[20], w[1], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[20], w[1], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[20], w[1], cxtype(cIPC[0], cIPC[1]), 0., 0., w[12]); 
  // Amplitude(s) for diagram number 75
  VVV1_0(w[7], w[4], w[12], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 76
  VVV1_0(w[1], w[7], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 77
  FFV1_0(w[3], w[11], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += +amp[0]; 
  jamp[14] += -amp[0]; 
  jamp[15] += +amp[0]; 
  jamp[16] += -amp[0]; 
  // Amplitude(s) for diagram number 78
  FFV1_0(w[21], w[11], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 79
  FFV1_0(w[15], w[2], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[10] += +amp[0]; 
  jamp[18] += -amp[0]; 
  // Amplitude(s) for diagram number 80
  FFV1_0(w[15], w[23], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[9], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[23]); 
  // Amplitude(s) for diagram number 81
  FFV1_0(w[15], w[23], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += -amp[0]; 
  FFV1_2(w[15], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  // Amplitude(s) for diagram number 82
  FFV1_0(w[12], w[9], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += -amp[0]; 
  // Amplitude(s) for diagram number 83
  FFV1_0(w[13], w[23], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += -amp[0]; 
  FFV1_2(w[13], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[21]); 
  // Amplitude(s) for diagram number 84
  FFV1_0(w[21], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += -amp[0]; 
  // Amplitude(s) for diagram number 85
  FFV1_0(w[3], w[23], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[10], cxtype(cIPC[0], cIPC[1]), 0., 0., w[23]); 
  // Amplitude(s) for diagram number 86
  FFV1_0(w[3], w[9], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[11] += +amp[0]; 
  FFV1_2(w[16], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[22]); 
  // Amplitude(s) for diagram number 87
  FFV1_0(w[22], w[11], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += -amp[0]; 
  FFV1_1(w[11], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[20]); 
  // Amplitude(s) for diagram number 88
  FFV1_0(w[16], w[20], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += -amp[0]; 
  // Amplitude(s) for diagram number 89
  FFV1_0(w[22], w[14], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[22] += -amp[0]; 
  FFV1_1(w[14], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[24]); 
  // Amplitude(s) for diagram number 90
  FFV1_0(w[16], w[24], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += -amp[0]; 
  // Amplitude(s) for diagram number 91
  FFV1_0(w[22], w[2], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 92
  FFV1_0(w[16], w[2], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +amp[0]; 
  jamp[5] += -amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[22] += +amp[0]; 
  // Amplitude(s) for diagram number 93
  VVVV1_0(w[0], w[6], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[0], w[6], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[0], w[6], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[6], cxtype(cIPC[0], cIPC[1]), 0., 0., w[22]); 
  // Amplitude(s) for diagram number 94
  VVV1_0(w[7], w[5], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[7], cxtype(cIPC[0], cIPC[1]), 0., 0., w[25]); 
  // Amplitude(s) for diagram number 95
  VVV1_0(w[6], w[5], w[25], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 96
  FFV1_0(w[3], w[14], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +amp[0]; 
  jamp[19] += -amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[23] += +amp[0]; 
  // Amplitude(s) for diagram number 97
  FFV1_0(w[3], w[24], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 98
  FFV1_0(w[13], w[2], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[14] += +amp[0]; 
  // Amplitude(s) for diagram number 99
  FFV1_0(w[21], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 100
  VVVV1_0(w[0], w[18], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[0], w[18], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[0], w[18], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[18], cxtype(cIPC[0], cIPC[1]), 0., 0., w[6]); 
  // Amplitude(s) for diagram number 101
  VVV1_0(w[7], w[4], w[6], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 102
  VVV1_0(w[18], w[4], w[25], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 103
  FFV1_0(w[3], w[11], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[17] += +amp[0]; 
  // Amplitude(s) for diagram number 104
  FFV1_0(w[3], w[20], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 105
  FFV1_0(w[15], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[4] += -amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[20] += +amp[0]; 
  // Amplitude(s) for diagram number 106
  FFV1_0(w[12], w[2], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 107
  VVVV1_0(w[0], w[1], w[7], w[10], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[0], w[1], w[7], w[10], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[0], w[1], w[7], w[10], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 108
  VVV1_0(w[1], w[10], w[25], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 109
  VVV1_0(w[1], w[7], w[23], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 110
  FFV1_0(w[13], w[20], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += -amp[0]; 
  // Amplitude(s) for diagram number 111
  FFV1_0(w[21], w[11], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[14] += -amp[0]; 
  // Amplitude(s) for diagram number 112
  FFV1_0(w[15], w[24], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += -amp[0]; 
  // Amplitude(s) for diagram number 113
  FFV1_0(w[12], w[14], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[20] += -amp[0]; 
  VVVV1P0_1(w[0], w[1], w[4], cxtype(cIPC[4], cIPC[5]), 0., 0., w[12]); 
  VVVV3P0_1(w[0], w[1], w[4], cxtype(cIPC[4], cIPC[5]), 0., 0., w[24]); 
  VVVV4P0_1(w[0], w[1], w[4], cxtype(cIPC[4], cIPC[5]), 0., 0., w[21]); 
  // Amplitude(s) for diagram number 114
  VVV1_0(w[12], w[7], w[5], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[24], w[7], w[5], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[21], w[7], w[5], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  jamp[23] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 115
  FFV1_0(w[3], w[14], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +amp[0]; 
  jamp[19] += -amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_0(w[3], w[14], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += -amp[0]; 
  jamp[20] += +amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[22] += +amp[0]; 
  FFV1_0(w[3], w[14], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += -amp[0]; 
  jamp[20] += +amp[0]; 
  jamp[22] += +amp[0]; 
  jamp[23] += -amp[0]; 
  // Amplitude(s) for diagram number 116
  FFV1_0(w[13], w[2], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[14] += +amp[0]; 
  FFV1_0(w[13], w[2], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += -amp[0]; 
  jamp[6] += +amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[12] += +amp[0]; 
  FFV1_0(w[13], w[2], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += -amp[0]; 
  jamp[6] += +amp[0]; 
  jamp[12] += +amp[0]; 
  jamp[14] += -amp[0]; 
  VVVV1P0_1(w[0], w[1], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[21]); 
  VVVV3P0_1(w[0], w[1], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[13]); 
  VVVV4P0_1(w[0], w[1], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[24]); 
  // Amplitude(s) for diagram number 117
  VVV1_0(w[21], w[7], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[13], w[7], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[24], w[7], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 118
  FFV1_0(w[3], w[11], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[17] += +amp[0]; 
  FFV1_0(w[3], w[11], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += -amp[0]; 
  jamp[14] += +amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[16] += +amp[0]; 
  FFV1_0(w[3], w[11], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += -amp[0]; 
  jamp[14] += +amp[0]; 
  jamp[16] += +amp[0]; 
  jamp[17] += -amp[0]; 
  // Amplitude(s) for diagram number 119
  FFV1_0(w[15], w[2], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[4] += -amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[20] += +amp[0]; 
  FFV1_0(w[15], w[2], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += -amp[0]; 
  jamp[7] += +amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[18] += +amp[0]; 
  FFV1_0(w[15], w[2], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += -amp[0]; 
  jamp[7] += +amp[0]; 
  jamp[18] += +amp[0]; 
  jamp[20] += -amp[0]; 
  VVVV1P0_1(w[0], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[24]); 
  VVVV3P0_1(w[0], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[15]); 
  VVVV4P0_1(w[0], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[13]); 
  // Amplitude(s) for diagram number 120
  FFV1_0(w[3], w[9], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[11] += +amp[0]; 
  FFV1_0(w[3], w[9], w[15], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += -amp[0]; 
  jamp[8] += +amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[10] += +amp[0]; 
  FFV1_0(w[3], w[9], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += -amp[0]; 
  jamp[8] += +amp[0]; 
  jamp[10] += +amp[0]; 
  jamp[11] += -amp[0]; 
  // Amplitude(s) for diagram number 121
  FFV1_0(w[16], w[2], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +amp[0]; 
  jamp[5] += -amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[22] += +amp[0]; 
  FFV1_0(w[16], w[2], w[15], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += -amp[0]; 
  jamp[13] += +amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[19] += +amp[0]; 
  FFV1_0(w[16], w[2], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += -amp[0]; 
  jamp[13] += +amp[0]; 
  jamp[19] += +amp[0]; 
  jamp[22] += -amp[0]; 
  // Amplitude(s) for diagram number 122
  VVV1_0(w[24], w[1], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[15], w[1], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[13], w[1], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 123
  VVV1_0(w[0], w[17], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += -cxtype(0, 1) * amp[0]; 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[23] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[0], w[19], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[0], w[8], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // double CPPProcess::matrix_1_gg_ttxgg() {



  // Local variables

  // The color matrix;
  static const double denom[ncolor] = {54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54};
  static const double cf[ncolor][ncolor] = {{512, -64, -64, 8, 8, 80, -64, 8,
      8, -1, -1, -10, 8, -1, 80, -10, 71, 62, -1, -10, -10, 62, 62, -28}, {-64,
      512, 8, 80, -64, 8, 8, -64, -1, -10, 8, -1, -1, -10, -10, 62, 62, -28, 8,
      -1, 80, -10, 71, 62}, {-64, 8, 512, -64, 80, 8, 8, -1, 80, -10, 71, 62,
      -64, 8, 8, -1, -1, -10, -10, -1, 62, -28, -10, 62}, {8, 80, -64, 512, 8,
      -64, -1, -10, -10, 62, 62, -28, 8, -64, -1, -10, 8, -1, -1, 8, 71, 62,
      80, -10}, {8, -64, 80, 8, 512, -64, -1, 8, 71, 62, 80, -10, -10, -1, 62,
      -28, -10, 62, -64, 8, 8, -1, -1, -10}, {80, 8, 8, -64, -64, 512, -10, -1,
      62, -28, -10, 62, -1, 8, 71, 62, 80, -10, 8, -64, -1, -10, 8, -1}, {-64,
      8, 8, -1, -1, -10, 512, -64, -64, 8, 8, 80, 80, -10, 8, -1, 62, 71, -10,
      62, -1, -10, -28, 62}, {8, -64, -1, -10, 8, -1, -64, 512, 8, 80, -64, 8,
      -10, 62, -1, -10, -28, 62, 80, -10, 8, -1, 62, 71}, {8, -1, 80, -10, 71,
      62, -64, 8, 512, -64, 80, 8, 8, -1, -64, 8, -10, -1, 62, -28, -10, -1,
      62, -10}, {-1, -10, -10, 62, 62, -28, 8, 80, -64, 512, 8, -64, -1, -10,
      8, -64, -1, 8, 71, 62, -1, 8, -10, 80}, {-1, 8, 71, 62, 80, -10, 8, -64,
      80, 8, 512, -64, 62, -28, -10, -1, 62, -10, 8, -1, -64, 8, -10, -1},
      {-10, -1, 62, -28, -10, 62, 80, 8, 8, -64, -64, 512, 71, 62, -1, 8, -10,
      80, -1, -10, 8, -64, -1, 8}, {8, -1, -64, 8, -10, -1, 80, -10, 8, -1, 62,
      71, 512, -64, -64, 8, 8, 80, 62, -10, -28, 62, -1, -10}, {-1, -10, 8,
      -64, -1, 8, -10, 62, -1, -10, -28, 62, -64, 512, 8, 80, -64, 8, -10, 80,
      62, 71, 8, -1}, {80, -10, 8, -1, 62, 71, 8, -1, -64, 8, -10, -1, -64, 8,
      512, -64, 80, 8, -28, 62, 62, -10, -10, -1}, {-10, 62, -1, -10, -28, 62,
      -1, -10, 8, -64, -1, 8, 8, 80, -64, 512, 8, -64, 62, 71, -10, 80, -1, 8},
      {71, 62, -1, 8, -10, 80, 62, -28, -10, -1, 62, -10, 8, -64, 80, 8, 512,
      -64, -1, 8, -10, -1, -64, 8}, {62, -28, -10, -1, 62, -10, 71, 62, -1, 8,
      -10, 80, 80, 8, 8, -64, -64, 512, -10, -1, -1, 8, 8, -64}, {-1, 8, -10,
      -1, -64, 8, -10, 80, 62, 71, 8, -1, 62, -10, -28, 62, -1, -10, 512, -64,
      -64, 8, 8, 80}, {-10, -1, -1, 8, 8, -64, 62, -10, -28, 62, -1, -10, -10,
      80, 62, 71, 8, -1, -64, 512, 8, 80, -64, 8}, {-10, 80, 62, 71, 8, -1, -1,
      8, -10, -1, -64, 8, -28, 62, 62, -10, -10, -1, -64, 8, 512, -64, 80, 8},
      {62, -10, -28, 62, -1, -10, -10, -1, -1, 8, 8, -64, 62, 71, -10, 80, -1,
      8, 8, 80, -64, 512, 8, -64}, {62, 71, -10, 80, -1, 8, -28, 62, 62, -10,
      -10, -1, -1, 8, -10, -1, -64, 8, 8, -64, 80, 8, 512, -64}, {-28, 62, 62,
      -10, -10, -1, 62, 71, -10, 80, -1, 8, -10, -1, -1, 8, 8, -64, 80, 8, 8,
      -64, -64, 512}};


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

  mgDebug(1, __FUNCTION__); 
  return; 
}



CPPProcess::CPPProcess(int numiterations, int gpublocks, int gputhreads, 
bool verbose, bool debug)
: m_numiterations(numiterations), gpu_nblocks(gpublocks), 
gpu_nthreads(gputhreads), m_verbose(verbose), 
dim(gpu_nblocks * gpu_nthreads) 
{


  // Helicities for the process - nodim
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1}, {-1, -1,
      -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1}, {-1, -1,
      -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1}, {-1, -1,
      -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1}, {-1, -1,
      1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1}, {-1, -1, 1,
      1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1, 1, -1, -1,
      -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1, 1, -1, -1,
      1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, 1,
      -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1, -1, -1, 1},
      {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1, -1, -1}, {-1,
      1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1}, {1, -1, -1,
      -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1}, {1, -1, -1,
      -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1, -1, -1, 1,
      1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, -1,
      1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1, -1, -1},
      {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1}, {1, 1,
      -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1, 1, -1,
      -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1, 1, 1,
      -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1, 1},
      {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1, 1,
      1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
#ifdef __CUDACC__
  checkCuda(cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int))); 
#else
  memcpy(cHel, tHel, ncomb * nexternal * sizeof(int)); 
#endif

  // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
  assert(sizeof(cxtype) == 2 * sizeof(fptype)); 
}

CPPProcess::~CPPProcess() {}

const std::vector<fptype> &CPPProcess::getMasses() const {return mME;}

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name) 
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
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 

  static cxtype tIPC[3] = {pars->GC_10, pars->GC_11, pars->GC_12}; 
  static double tIPD[2] = {pars->mdl_MT, pars->mdl_WT}; 

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
  //** START LOOP ON IEVT **
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
    const int denominators[1] = {512}; 

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
    //** END LOOP ON IEVT **
    mgDebugFinalise(); 

  }

  //--------------------------------------------------------------------------

}



}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


