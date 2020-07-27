//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

//#include "HelAmps_sm.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
#include <thrust/complex.h> 
using namespace std; 

namespace MG5_sm 
{

__device__ 
void ixxxxx( const double pvec[3], double fmass, int nhel, int nsf, thrust::complex<double> fi[6] ) 
{
  thrust::complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 

  double p[4] = {0, pvec[0], pvec[1], pvec[2]}; 
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + fmass * fmass); 
  fi[0] = thrust::complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = thrust::complex<double> (-p[1] * nsf, -p[2] * nsf); 
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
      chi[0] = thrust::complex<double> (sqrt(pp3 * 0.5/pp), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = thrust::complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = 
        thrust::complex<double> (nh * p[1], p[2])/sqrt(2.0 * pp * pp3); 
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
    chi[0] = thrust::complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = thrust::complex<double> (-nhel * sqrt(2.0 * p[0]), 0.0); 
    }
    else
    {
      chi[1] = thrust::complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = thrust::complex<double> (0.0, 0.0); 
      fi[3] = thrust::complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = thrust::complex<double> (0.0, 0.0); 
      fi[5] = thrust::complex<double> (0.0, 0.0); 
    }
  }

  return; 
}

__device__ 
void txxxxx( const double pvec[3], double tmass, int nhel, int nst, thrust::complex<double> tc[18] ) 
{
  thrust::complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  double p[4] = {0, pvec[0], pvec[1], pvec[2]}; 
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + tmass * tmass); 
  sqh = sqrt(0.5); 
  sqs = sqrt(0.5/3); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], sqrt(pt2 + p[3] * p[3])); 
  pt = min(pp, sqrt(pt2)); 

  ft[4][0] = thrust::complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = thrust::complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if (nhel >= 0)
  {
    if (pp == 0)
    {
      ep[0] = thrust::complex<double> (0, 0); 
      ep[1] = thrust::complex<double> (-sqh, 0); 
      ep[2] = thrust::complex<double> (0, nst * sqh); 
      ep[3] = thrust::complex<double> (0, 0); 
    }
    else
    {
      ep[0] = thrust::complex<double> (0, 0); 
      ep[3] = thrust::complex<double> (pt/pp * sqh, 0); 

      if (pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = thrust::complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = thrust::complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = thrust::complex<double> (-sqh, 0); 
        ep[2] = 
        thrust::complex<double> (0, nst * (p[3] < 0) ? - abs(sqh) : abs(sqh)); 
      }
    }
  }

  // construct eps-
  if (nhel <= 0)
  {
    if (pp == 0)
    {
      em[0] = thrust::complex<double> (0, 0); 
      em[1] = thrust::complex<double> (sqh, 0); 
      em[2] = thrust::complex<double> (0, nst * sqh); 
      em[3] = thrust::complex<double> (0, 0); 
    }
    else
    {
      em[0] = thrust::complex<double> (0, 0); 
      em[3] = thrust::complex<double> (-pt/pp * sqh, 0); 

      if (pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = thrust::complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = thrust::complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = thrust::complex<double> (sqh, 0); 
        em[2] = 
        thrust::complex<double> (0, nst * (p[3] < 0) ? - abs(sqh) : abs(sqh)); 
      }
    }
  }

  // construct eps0
  if (std::labs(nhel) <= 1)
  {
    if (pp == 0)
    {
      e0[0] = thrust::complex<double> (0, 0); 
      e0[1] = thrust::complex<double> (0, 0); 
      e0[2] = thrust::complex<double> (0, 0); 
      e0[3] = thrust::complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = thrust::complex<double> (pp/tmass, 0); 
      e0[3] = thrust::complex<double> (p[3] * emp, 0); 

      if (pt != 0)
      {
        e0[1] = thrust::complex<double> (p[1] * emp, 0); 
        e0[2] = thrust::complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = thrust::complex<double> (0, 0); 
        e0[2] = thrust::complex<double> (0, 0); 
      }
    }
  }

  if (nhel == 2)
  {
    for (j = 0; j < 4; j++ )
    {
      for (i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if (nhel == -2)
  {
    for (j = 0; j < 4; j++ )
    {
      for (i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if (tmass == 0)
  {
    for (j = 0; j < 4; j++ )
    {
      for (i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if (tmass != 0)
  {
    if (nhel == 1)
    {
      for (j = 0; j < 4; j++ )
      {
        for (i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if (nhel == 0)
    {
      for (j = 0; j < 4; j++ )
      {
        for (i = 0; i < 4; i++ )
          ft[i][j] = 
        sqs * (ep[i] * em[j] + em[i] * ep[j] + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if (nhel == -1)
    {
      for (j = 0; j < 4; j++ )
      {
        for (i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      // sr fixme // std::cerr << "Invalid helicity in txxxxx.\n";
      // sr fixme // std::exit(1);
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for (j = 0; j < 4; j++ )
  {
    for (i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

__device__ 
void vxxxxx( const double pvec[3], double vmass, int nhel, int nsv, thrust::complex<double> vc[6] ) 
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 

  double p[4] = {0, pvec[0], pvec[1], pvec[2]}; 
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + vmass * vmass); 

  sqh = sqrt(0.5); 
  hel = double(nhel); 
  nsvahl = nsv * std::abs(hel); 
  pt2 = (p[1] * p[1]) + (p[2] * p[2]); 
  pp = min(p[0], sqrt(pt2 + (p[3] * p[3]))); 
  pt = min(pp, sqrt(pt2)); 
  vc[0] = thrust::complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = thrust::complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - std::abs(hel); 
    if (pp == 0.0)
    {
      vc[2] = thrust::complex<double> (0.0, 0.0); 
      vc[3] = thrust::complex<double> (-hel * sqh, 0.0); 
      vc[4] = thrust::complex<double> (0.0, nsvahl * sqh); 
      vc[5] = thrust::complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = thrust::complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = 
      thrust::complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = thrust::complex<double> (hel0 * p[1] * emp - p[1] * pzpt, 
         - nsvahl * p[2]/pt * sqh); 
        vc[4] = thrust::complex<double> (hel0 * p[2] * emp - p[2] * pzpt, 
        nsvahl * p[1]/pt * sqh); 
      }
      else
      {
        vc[3] = thrust::complex<double> (-hel * sqh, 0.0); 
        vc[4] = thrust::complex<double> (0.0, nsvahl * (p[3] < 0) ? - abs(sqh)
        : abs(sqh)); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = sqrt((p[1] * p[1]) + (p[2] * p[2])); 
    vc[2] = thrust::complex<double> (0.0, 0.0); 
    vc[5] = thrust::complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = thrust::complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = thrust::complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = thrust::complex<double> (-hel * sqh, 0.0); 
      vc[4] = 
      thrust::complex<double> (0.0, nsv * (p[3] < 0) ? - abs(sqh) : abs(sqh)); 
    }
  }
  return; 
}

__device__ 
void sxxxxx( const double pvec[3], int nss, thrust::complex<double> sc[3] ) 
{
  // double p[4] = {0, pvec[0], pvec[1], pvec[2]};
  // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+fmass*fmass);
  double p[4] = {0, 0, 0, 0}; 
  printf("scalar not supported so far. to do: fix mass issue"); 
  sc[2] = thrust::complex<double> (1.00, 0.00); 
  sc[0] = thrust::complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = thrust::complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

__device__ 
void oxxxxx( const double pvec[3], double fmass, int nhel, int nsf, thrust::complex<double> fo[6] ) 
{
  thrust::complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 

  double p[4] = {0, pvec[0], pvec[1], pvec[2]}; 
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + fmass * fmass); 

  fo[0] = thrust::complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = thrust::complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
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
      pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = thrust::complex<double> (sqrt(pp3 * 0.5/pp), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = thrust::complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = 
        thrust::complex<double> (nh * p[1], -p[2])/sqrt(2.0 * pp * pp3); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if ((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.00)) * nsf; 
    }
    chi[0] = thrust::complex<double> (sqp0p3, 0.00); 
    if (sqp0p3 == 0.000)
    {
      chi[1] = thrust::complex<double> (-nhel, 0.00) * sqrt(2.0 * p[0]); 
    }
    else
    {
      chi[1] = thrust::complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = thrust::complex<double> (0.00, 0.00); 
      fo[5] = thrust::complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = thrust::complex<double> (0.00, 0.00); 
      fo[3] = thrust::complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}


__device__ 
void FFV2_0( const thrust::complex<double> F1[], 
             const thrust::complex<double> F2[], 
             const thrust::complex<double> V3[], 
             const thrust::complex<double> COUP, 
             thrust::complex<double> * vertex )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  thrust::complex<double> TMP0; 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))); 
  (*vertex) = COUP * - cI * TMP0; 
}


__device__ 
void FFV2_3( const thrust::complex<double> F1[], 
             const thrust::complex<double> F2[], 
             const thrust::complex<double> COUP, 
             const double M3, 
             const double W3, 
             thrust::complex<double> V3[] )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  thrust::complex<double> TMP1; 
  thrust::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) + 
  F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3]))); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * 
  P3[3]) - M3 * (M3 - cI * W3)); 
  V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP1); 
  V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * 
  TMP1); 
  V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2]
   * OM3 * TMP1); 
  V3[5] = denom * (-cI) * (-F1[2] * F2[4] - P3[3] * OM3 * TMP1 + F1[3] * 
  F2[5]); 
}


__device__ 
void FFV4_0( const thrust::complex<double> F1[], 
             const thrust::complex<double> F2[], 
             const thrust::complex<double> V3[], 
             const thrust::complex<double> COUP, 
             thrust::complex<double> * vertex )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  thrust::complex<double> TMP0; 
  thrust::complex<double> TMP2; 
  TMP2 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
  F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))); 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))); 
  (*vertex) = COUP * (-1.) * (+cI * (TMP0) + 2. * cI * (TMP2)); 
}


__device__ 
void FFV4_3( const thrust::complex<double> F1[], 
             const thrust::complex<double> F2[], 
             const thrust::complex<double> COUP, 
             const double M3, 
             const double W3, 
             thrust::complex<double> V3[] )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  thrust::complex<double> TMP1; 
  thrust::complex<double> TMP3; 
  thrust::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) + 
  F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3]))); 
  TMP3 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) + 
  F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3]))); 
  denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * 
  P3[3]) - M3 * (M3 - cI * W3)); 
  V3[2] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[0] * (TMP1 + 2. * (TMP3)) + 
  (+1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] * 
  F2[3])); 
  V3[3] = denom * (-2. * cI) * (OM3 * - 1./2. * P3[1] * (TMP1 + 2. * (TMP3)) + 
  (-1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] * 
  F2[2])); 
  V3[4] = denom * 2. * cI * (OM3 * 1./2. * P3[2] * (TMP1 + 2. * (TMP3)) + 
  (+1./2. * cI * (F1[2] * F2[5]) - 1./2. * cI * (F1[3] * F2[4]) - cI * 
  (F1[4] * F2[3]) + cI * (F1[5] * F2[2]))); 
  V3[5] = denom * 2. * cI * (OM3 * 1./2. * P3[3] * (TMP1 + 2. * (TMP3)) + 
  (+1./2. * (F1[2] * F2[4]) - 1./2. * (F1[3] * F2[5]) - F1[4] * F2[2] + 
  F1[5] * F2[3])); 
}


__device__ 
void FFV1_0( const thrust::complex<double> F1[], 
             const thrust::complex<double> F2[], 
             const thrust::complex<double> V3[], 
             const thrust::complex<double> COUP, 
             thrust::complex<double> * vertex )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  thrust::complex<double> TMP4; 
  TMP4 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) + 
  (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
  F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))))); 
  (*vertex) = COUP * - cI * TMP4; 
}


__device__ 
void FFV1P0_3( const thrust::complex<double> F1[], 
               const thrust::complex<double> F2[], 
               const thrust::complex<double> COUP, 
               const double M3, 
               const double W3, 
               thrust::complex<double> V3[])
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  double P3[4]; 
  thrust::complex<double> denom; 
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


__device__ 
void FFV2_4_0( const thrust::complex<double> F1[], 
               const thrust::complex<double> F2[], 
               const thrust::complex<double> V3[], 
               const thrust::complex<double> COUP1, 
               const thrust::complex<double> COUP2, 
               thrust::complex<double> * vertex )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  thrust::complex<double> TMP0; 
  thrust::complex<double> TMP2; 
  TMP2 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
  F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))); 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))); 
  (*vertex) = (-1.) * (COUP2 * (+cI * (TMP0) + 2. * cI * (TMP2)) + cI * (TMP0 * 
  COUP1)); 
}


__device__ 
void FFV2_4_3( const thrust::complex<double> F1[], 
               const thrust::complex<double> F2[], 
               const thrust::complex<double> COUP1, 
               const thrust::complex<double> COUP2, 
               const double M3, 
               const double W3, 
               thrust::complex<double> V3[] )
{
  thrust::complex<double> cI = thrust::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  thrust::complex<double> TMP1; 
  thrust::complex<double> TMP3; 
  thrust::complex<double> denom; 
  OM3 = 0.; 
  if (M3 != 0.)
    OM3 = 1./(M3 * M3); 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) + 
  F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3]))); 
  TMP3 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) + 
  F1[5] * (F2[2] * (-P3[1] + cI * (P3[2])) + F2[3] * (P3[0] + P3[3]))); 
  denom = 1./((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * 
  P3[3]) - M3 * (M3 - cI * W3)); 
  V3[2] = denom * (-2. * cI) * (COUP2 * (OM3 * - 1./2. * P3[0] * (TMP1 + 2. * 
  (TMP3)) + (+1./2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + 
  F1[5] * F2[3])) + 1./2. * (COUP1 * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0]
   * OM3 * TMP1))); 
  V3[3] = denom * (-2. * cI) * (COUP2 * (OM3 * - 1./2. * P3[1] * (TMP1 + 2. * 
  (TMP3)) + (-1./2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + 
  F1[5] * F2[2])) - 1./2. * (COUP1 * (F1[2] * F2[5] + F1[3] * F2[4] + P3[1]
   * OM3 * TMP1))); 
  V3[4] = denom * cI * (COUP2 * (OM3 * P3[2] * (TMP1 + 2. * (TMP3)) + (+cI * 
  (F1[2] * F2[5]) - cI * (F1[3] * F2[4]) - 2. * cI * (F1[4] * F2[3]) + 2. * 
  cI * (F1[5] * F2[2]))) + COUP1 * (+cI * (F1[2] * F2[5]) - cI * (F1[3] * 
  F2[4]) + P3[2] * OM3 * TMP1)); 
  V3[5] = denom * 2. * cI * (COUP2 * (OM3 * 1./2. * P3[3] * (TMP1 + 2. * 
  (TMP3)) + (+1./2. * (F1[2] * F2[4]) - 1./2. * (F1[3] * F2[5]) - F1[4] * 
  F2[2] + F1[5] * F2[3])) + 1./2. * (COUP1 * (F1[2] * F2[4] + P3[3] * OM3 * 
  TMP1 - F1[3] * F2[5]))); 
}


}  // end namespace $(namespace)s_sm


//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
//#include "HelAmps_sm.h"

#include <algorithm> 
#include <iostream> 
#include <thrust/complex.h> 

using namespace MG5_sm; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

__constant__ int cHel[16][4]; 
// __constant__ double cmME[4]; value hardcoded now
// extern __constant__ int cPerm[4];
// 
__constant__ double cIPC[6];  // coupling ?
__constant__ double cIPD[2]; 


__device__ void printComplex6( const thrust::complex<double> c6[6], const char* txt )
{
  for (int i=0; i<6; i++)
    printf("%s[%i]=(%f,%f)\n", txt, i, c6[i].real(), c6[i].imag());
}

__device__ void printDouble4( const double p[4], const char* txt )
{
  //for (int i=0; i<4; i++) printf("%s[%i]=%f\n", txt, i, p[i]);
  printf("%s[0..3] = ", txt);
  for (int i=0; i<4; i++) printf("%f ", p[i]);
  printf("\n");
}

// Evaluate |M|^2 for each subprocess

__device__
void calculate_wavefunctions( int ihel, const double local_mom[4][3], double &matrix )
{
  thrust::complex<double> amp[2]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[5][6]; 
  oxxxxx(local_mom[0], 0., cHel[ihel][0], -1, w[0]); 
  ixxxxx(local_mom[1], 0., cHel[ihel][1], +1, w[1]); 
  ixxxxx(local_mom[2], 0., cHel[ihel][2], -1, w[2]); 
  oxxxxx(local_mom[3], 0., cHel[ihel][3], +1, w[3]); 
  FFV1P0_3(w[1], w[0], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[4]);
  // Amplitude(s) for diagram number 1
  FFV1_0(w[2], w[3], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[0]);
  FFV2_4_3(w[1], w[0], thrust::complex<double> (cIPC[2], cIPC[3]),
      thrust::complex<double> (cIPC[4], cIPC[5]), cIPD[0], cIPD[1], w[4]);
  // Amplitude(s) for diagram number 2
  FFV2_4_0(w[2], w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      thrust::complex<double> (cIPC[4], cIPC[5]), &amp[1]);

  //printf("DEBUG0: amp0=(%f,%f) amp1=(%f,%f)\n", amp[0].real(), amp[0].imag(), amp[1].real(), amp[1].imag());

  // double CPPProcess::matrix_1_epem_mupmum() {
  int i, j; 
  // Local variables

  // const int ngraphs = 2;
  const int ncolor = 1; 
  thrust::complex<double> ztemp; 
  thrust::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{1}}; 

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1]; 

  // Sum and square the color flows to get the matrix element
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + (ztemp * conj(jamp[i])).real()/denom[i]; 
  }

  //printf("DEBUG0: ME=%f\n", matrix);

  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

}



CPPProcess::CPPProcess( int numiterations, int gpublocks, int gputhreads, bool verbose, bool debug )
  : m_numiterations(numiterations)
  , gpu_nblocks(gpublocks)
  , gpu_nthreads(gputhreads)
  , m_verbose(verbose)
  , m_debug(debug)
  , dim(gpu_nblocks * gpu_nthreads) 
{


  // Helicities for the process - nodim
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1, -1, 1},
      {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1}, {-1, 1,
      1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1, 1, -1},
      {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1, 1, 1,
      1}};
  cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int)); 
  // perm - nodim
  // static int perm[nexternal] = {0, 1, 2, 3};
}

CPPProcess::~CPPProcess() {}

const std::vector<double> &CPPProcess::getMasses() const {return mME;}

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name) 
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
  static thrust::complex<double> tIPC[3] = {pars->GC_3, pars->GC_50,
      pars->GC_59};
  static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ}; 

  cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(thrust::complex<double> )); 
  cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(double)); 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

__global__
void sigmaKin( const double * allmomenta, double * output ) 
{
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();

  // Reset color flows

  // for (int xx = 0; xx < 384; ++xx) {
  const int nprocesses = 1; 
  int tid = blockIdx.x * blockDim.x + threadIdx.x; 

  // char *devPtr = (char *)tp.ptr;
  // size_t dpt = tp.pitch;
  // size_t slicePitch = dpt * 4;

  // char *dps = devPtr + dim * slicePitch;
  double matrix_element[nprocesses]; 

  thrust::complex<double> amp[2]; 

  double local_m[4][3]; 
  int DIM = blockDim.x * gridDim.x; 
  // for (int i=0; i<20;i++){
  // printf(" %f ", allmomenta[i]);
  // }
  // printf("\n");
  // printf("DIM is %i/%i\n", tid, DIM);
  for (int i = 0; i < 4; i++ )
  {
    for (int j = 0; j < 3; j++ )
    {
      local_m[i][j] = allmomenta[i * 3 * DIM + j * DIM + tid]; 
      // printf(" %f ", local_m[i][j]);
    }
    // printf("\n");
  }


  // Local variables and constants
  const int ncomb = 16; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  // double t[1];
  // Helicities for the process
  // static const int helicities[ncomb][nexternal] =
  // {{-1,-1,-1,-1},{-1,-1,-1,1},{-1,-1,1,-1},{-1,-1,1,1},{-1,1,-1,-1},{-1,1,-1,
  // 1},{-1,1,1,-1},{-1,1,1,1},{1,-1,-1,-1},{1,-1,-1,1},{1,-1,1,-1},{1,-1,1,1},{
  // 1,1,-1,-1},{1,1,-1,1},{1,1,1,-1},{1,1,1,1}};
  // Denominators: spins, colors and identical particles
  const int denominators[1] = {4}; 


  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  // int perm[nexternal];
  // for(int i = 0; i < nexternal; i++){
  // perm[i]=i;
  // }


  //printf("\n\nblock (%i / %i), thread (%i)\n\n", blockIdx.x, blockDim.x, threadIdx.x);
  for (int ihel = 0; ihel < ncomb; ihel++ )
  {
    calculate_wavefunctions(ihel, local_m, matrix_element[0]); 
  }


  for (int i = 0; i < nprocesses; ++ i)
  {
    matrix_element[i] /= denominators[i]; 
  }
  for (int i = 0; i < nprocesses; ++ i)
  {
    output[i * nprocesses + tid] = matrix_element[i]; 
    // printf("output %i %i %i %f", tid, i, i*nprocesses+tid,
    // output[i*nprocesses+tid]);

  }
  //printf("DEBUG: ME=%f\n", matrix_element[0]);

}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


