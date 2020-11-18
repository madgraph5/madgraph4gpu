//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H


//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================




#include "Parameters_sm.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Complex.hpp"

template<typename Stride, typename Device, typename MemSpace>
KOKKOS_FUNCTION void ixxxxx(Kokkos::View<double*,Stride,Device,MemSpace> pvec, double fmass, int nhel, int nsf, 
Kokkos::complex<double> fi[6])
{
  Kokkos::complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 

  double p[4] = {0, pvec(1), pvec(2), pvec(3)};
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + fmass * fmass); 
  fi[0] = Kokkos::complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = Kokkos::complex<double> (-p[1] * nsf, -p[2] * nsf); 
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
      chi[0] = Kokkos::complex<double> (sqrt(pp3 * 0.5/pp), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = Kokkos::complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = 
        Kokkos::complex<double> (nh * p[1], p[2])/sqrt(2.0 * pp * pp3); 
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
    chi[0] = Kokkos::complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = Kokkos::complex<double> (-nhel * sqrt(2.0 * p[0]), 0.0); 
    }
    else
    {
      chi[1] = Kokkos::complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = Kokkos::complex<double> (0.0, 0.0); 
      fi[3] = Kokkos::complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = Kokkos::complex<double> (0.0, 0.0); 
      fi[5] = Kokkos::complex<double> (0.0, 0.0); 
    }
  }

  return; 
}


template<typename Stride, typename Device, typename MemSpace>
KOKKOS_FUNCTION void txxxxx(Kokkos::View<double*,Stride,Device,MemSpace> pvec, double tmass, int nhel, int nst, 
Kokkos::complex<double> tc[18]) 
{
  Kokkos::complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  double p[4] = {0, pvec(1), pvec(2), pvec(3)}; 
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + tmass * tmass); 
  sqh = sqrt(0.5); 
  sqs = sqrt(0.5/3); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], sqrt(pt2 + p[3] * p[3])); 
  pt = min(pp, sqrt(pt2)); 

  ft[4][0] = Kokkos::complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = Kokkos::complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if (nhel >= 0)
  {
    if (pp == 0)
    {
      ep[0] = Kokkos::complex<double> (0, 0); 
      ep[1] = Kokkos::complex<double> (-sqh, 0); 
      ep[2] = Kokkos::complex<double> (0, nst * sqh); 
      ep[3] = Kokkos::complex<double> (0, 0); 
    }
    else
    {
      ep[0] = Kokkos::complex<double> (0, 0); 
      ep[3] = Kokkos::complex<double> (pt/pp * sqh, 0); 

      if (pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = Kokkos::complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = Kokkos::complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = Kokkos::complex<double> (-sqh, 0); 
        ep[2] = 
        Kokkos::complex<double> (0, nst * (p[3] < 0) ? - abs(sqh) : abs(sqh)); 
      }
    }
  }

  // construct eps-
  if (nhel <= 0)
  {
    if (pp == 0)
    {
      em[0] = Kokkos::complex<double> (0, 0); 
      em[1] = Kokkos::complex<double> (sqh, 0); 
      em[2] = Kokkos::complex<double> (0, nst * sqh); 
      em[3] = Kokkos::complex<double> (0, 0); 
    }
    else
    {
      em[0] = Kokkos::complex<double> (0, 0); 
      em[3] = Kokkos::complex<double> (-pt/pp * sqh, 0); 

      if (pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = Kokkos::complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = Kokkos::complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = Kokkos::complex<double> (sqh, 0); 
        em[2] = 
        Kokkos::complex<double> (0, nst * (p[3] < 0) ? - abs(sqh) : abs(sqh)); 
      }
    }
  }

  // construct eps0
  if (std::labs(nhel) <= 1)
  {
    if (pp == 0)
    {
      e0[0] = Kokkos::complex<double> (0, 0); 
      e0[1] = Kokkos::complex<double> (0, 0); 
      e0[2] = Kokkos::complex<double> (0, 0); 
      e0[3] = Kokkos::complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = Kokkos::complex<double> (pp/tmass, 0); 
      e0[3] = Kokkos::complex<double> (p[3] * emp, 0); 

      if (pt != 0)
      {
        e0[1] = Kokkos::complex<double> (p[1] * emp, 0); 
        e0[2] = Kokkos::complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = Kokkos::complex<double> (0, 0); 
        e0[2] = Kokkos::complex<double> (0, 0); 
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


template<typename Stride, typename Device, typename MemSpace>
KOKKOS_FUNCTION void vxxxxx(Kokkos::View<double*,Stride,Device,MemSpace> pvec, 
    double vmass, int nhel, int nsv, Kokkos::complex<double> vc[6]) 
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 

  double p[4] = {0, pvec(1), pvec(2), pvec(3)}; 
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + vmass * vmass); 

  sqh = sqrt(0.5); 
  hel = double(nhel); 
  nsvahl = nsv * std::abs(hel); 
  pt2 = (p[1] * p[1]) + (p[2] * p[2]); 
  pp = min(p[0], sqrt(pt2 + (p[3] * p[3]))); 
  pt = min(pp, sqrt(pt2)); 
  vc[0] = Kokkos::complex<double> (p[0] * nsv, p[3] * nsv);
  vc[1] = Kokkos::complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - std::abs(hel); 
    if (pp == 0.0)
    {
      vc[2] = Kokkos::complex<double> (0.0, 0.0); 
      vc[3] = Kokkos::complex<double> (-hel * sqh, 0.0); 
      vc[4] = Kokkos::complex<double> (0.0, nsvahl * sqh); 
      vc[5] = Kokkos::complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = Kokkos::complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = 
      Kokkos::complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = Kokkos::complex<double> (hel0 * p[1] * emp - p[1] * pzpt, 
         - nsvahl * p[2]/pt * sqh); 
        vc[4] = Kokkos::complex<double> (hel0 * p[2] * emp - p[2] * pzpt, 
        nsvahl * p[1]/pt * sqh); 
      }
      else
      {
        vc[3] = Kokkos::complex<double> (-hel * sqh, 0.0); 
        vc[4] = Kokkos::complex<double> (0.0, nsvahl * (p[3] < 0) ? - abs(sqh)
        : abs(sqh)); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = sqrt((p[1] * p[1]) + (p[2] * p[2])); 
    vc[2] = Kokkos::complex<double> (0.0, 0.0); 
    vc[5] = Kokkos::complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = Kokkos::complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = Kokkos::complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    { 
      vc[3] = Kokkos::complex<double> (-hel * sqh, 0.0);  
      vc[4] = 
      Kokkos::complex<double> (0.0, nsv * (p[3] < 0) ? - abs(sqh) : abs(sqh)); 
    }
  }
  return; 
}


template<typename Stride, typename Device, typename MemSpace>
KOKKOS_FUNCTION void sxxxxx(Kokkos::View<double*,Stride,Device,MemSpace> pvec, int nss, Kokkos::complex<double> sc[3]) 
{
  // double p[4] = {0, pvec[0], pvec[1], pvec[2]};
  // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+fmass*fmass);
  double p[4] = {0, 0, 0, 0}; 
  printf("scalar not supported so far. to do: fix mass issue"); 
  sc[2] = Kokkos::complex<double> (1.00, 0.00); 
  sc[0] = Kokkos::complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = Kokkos::complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}


template<typename Stride, typename Device, typename MemSpace>
KOKKOS_FUNCTION void oxxxxx(Kokkos::View<double*,Stride,Device,MemSpace> pvec, double fmass, int nhel, int nsf, 
Kokkos::complex<double> fo[6]) 
{
  Kokkos::complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 

  double p[4] = {0, pvec(1), pvec(2), pvec(3)};
  p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3] + fmass * fmass); 

  fo[0] = Kokkos::complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = Kokkos::complex<double> (p[1] * nsf, p[2] * nsf); 
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
      chi[0] = Kokkos::complex<double> (sqrt(pp3 * 0.5/pp), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = Kokkos::complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = 
        Kokkos::complex<double> (nh * p[1], -p[2])/sqrt(2.0 * pp * pp3); 
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
    chi[0] = Kokkos::complex<double> (sqp0p3, 0.00); 
    if (sqp0p3 == 0.000)
    {
      chi[1] = Kokkos::complex<double> (-nhel, 0.00) * sqrt(2.0 * p[0]); 
    }
    else
    {
      chi[1] = Kokkos::complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = Kokkos::complex<double> (0.00, 0.00); 
      fo[5] = Kokkos::complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = Kokkos::complex<double> (0.00, 0.00); 
      fo[3] = Kokkos::complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}


KOKKOS_FUNCTION void VVVV3_0(Kokkos::complex<double> V1[], const
Kokkos::complex<double> V2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> V4[], const Kokkos::complex<double> COUP, 
Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  Kokkos::complex<double> TMP0; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP2; 
  Kokkos::complex<double> TMP3; 
  TMP3 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP1 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP2 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP0 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  (*vertex) = COUP * (-cI * (TMP0 * TMP1) + cI * (TMP2 * TMP3)); 
}


KOKKOS_FUNCTION void VVVV3P0_1(Kokkos::complex<double> V2[], const
Kokkos::complex<double> V3[], const Kokkos::complex<double> V4[], const
Kokkos::complex<double> COUP, const double M1, const double W1, 
Kokkos::complex<double> V1[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P1[4]; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP2; 
  Kokkos::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP1 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP2 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (-cI * (V4[2] * TMP1) + cI * (V2[2] * TMP2)); 
  V1[3] = denom * (-cI * (V4[3] * TMP1) + cI * (V2[3] * TMP2)); 
  V1[4] = denom * (-cI * (V4[4] * TMP1) + cI * (V2[4] * TMP2)); 
  V1[5] = denom * (-cI * (V4[5] * TMP1) + cI * (V2[5] * TMP2)); 
}


KOKKOS_FUNCTION void VVVV1_0(Kokkos::complex<double> V1[], const
Kokkos::complex<double> V2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> V4[], const Kokkos::complex<double> COUP, 
Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  Kokkos::complex<double> TMP0; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP4; 
  Kokkos::complex<double> TMP5; 
  TMP5 = (V1[2] * V3[2] - V1[3] * V3[3] - V1[4] * V3[4] - V1[5] * V3[5]); 
  TMP1 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP0 = (V4[2] * V1[2] - V4[3] * V1[3] - V4[4] * V1[4] - V4[5] * V1[5]); 
  TMP4 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  (*vertex) = COUP * (-cI * (TMP0 * TMP1) + cI * (TMP4 * TMP5)); 
}


KOKKOS_FUNCTION void VVVV1P0_1(Kokkos::complex<double> V2[], const
Kokkos::complex<double> V3[], const Kokkos::complex<double> V4[], const
Kokkos::complex<double> COUP, const double M1, const double W1, 
Kokkos::complex<double> V1[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P1[4]; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP4; 
  Kokkos::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP1 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP4 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (-cI * (V4[2] * TMP1) + cI * (V3[2] * TMP4)); 
  V1[3] = denom * (-cI * (V4[3] * TMP1) + cI * (V3[3] * TMP4)); 
  V1[4] = denom * (-cI * (V4[4] * TMP1) + cI * (V3[4] * TMP4)); 
  V1[5] = denom * (-cI * (V4[5] * TMP1) + cI * (V3[5] * TMP4)); 
}


KOKKOS_FUNCTION void FFV1_0(Kokkos::complex<double> F1[], const
Kokkos::complex<double> F2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> COUP, Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  Kokkos::complex<double> TMP6; 
  TMP6 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) + 
  (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
  F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))))); 
  (*vertex) = COUP * - cI * TMP6; 
}


KOKKOS_FUNCTION void FFV1_1(Kokkos::complex<double> F2[], const
Kokkos::complex<double> V3[], const Kokkos::complex<double> COUP, const
double M1, const double W1, Kokkos::complex<double> F1[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P1[4]; 
  Kokkos::complex<double> denom; 
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


KOKKOS_FUNCTION void FFV1_2(Kokkos::complex<double> F1[], const
Kokkos::complex<double> V3[], const Kokkos::complex<double> COUP, const
double M2, const double W2, Kokkos::complex<double> F2[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P2[4]; 
  Kokkos::complex<double> denom; 
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


KOKKOS_FUNCTION void FFV1P0_3(Kokkos::complex<double> F1[], const
Kokkos::complex<double> F2[], const Kokkos::complex<double> COUP, const
double M3, const double W3, Kokkos::complex<double> V3[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P3[4]; 
  Kokkos::complex<double> denom; 
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


KOKKOS_FUNCTION void VVVV4_0(Kokkos::complex<double> V1[], const
Kokkos::complex<double> V2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> V4[], const Kokkos::complex<double> COUP, 
Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  Kokkos::complex<double> TMP2; 
  Kokkos::complex<double> TMP3; 
  Kokkos::complex<double> TMP4; 
  Kokkos::complex<double> TMP5; 
  TMP3 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP5 = (V1[2] * V3[2] - V1[3] * V3[3] - V1[4] * V3[4] - V1[5] * V3[5]); 
  TMP2 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP4 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  (*vertex) = COUP * (-cI * (TMP4 * TMP5) + cI * (TMP2 * TMP3)); 
}


KOKKOS_FUNCTION void VVVV4P0_1(Kokkos::complex<double> V2[], const
Kokkos::complex<double> V3[], const Kokkos::complex<double> V4[], const
Kokkos::complex<double> COUP, const double M1, const double W1, 
Kokkos::complex<double> V1[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P1[4]; 
  Kokkos::complex<double> TMP2; 
  Kokkos::complex<double> TMP4; 
  Kokkos::complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP2 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP4 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (-cI * (V3[2] * TMP4) + cI * (V2[2] * TMP2)); 
  V1[3] = denom * (-cI * (V3[3] * TMP4) + cI * (V2[3] * TMP2)); 
  V1[4] = denom * (-cI * (V3[4] * TMP4) + cI * (V2[4] * TMP2)); 
  V1[5] = denom * (-cI * (V3[5] * TMP4) + cI * (V2[5] * TMP2)); 
}


KOKKOS_FUNCTION void VVV1_0(Kokkos::complex<double> V1[], const
Kokkos::complex<double> V2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> COUP, Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP10; 
  Kokkos::complex<double> TMP11; 
  Kokkos::complex<double> TMP12; 
  Kokkos::complex<double> TMP3; 
  Kokkos::complex<double> TMP5; 
  Kokkos::complex<double> TMP7; 
  Kokkos::complex<double> TMP8; 
  Kokkos::complex<double> TMP9; 
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
  TMP9 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  TMP8 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP3 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  TMP1 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP7 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP5 = (V1[2] * V3[2] - V1[3] * V3[3] - V1[4] * V3[4] - V1[5] * V3[5]); 
  TMP10 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP11 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP12 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  (*vertex) = COUP * (TMP1 * (-cI * (TMP11) + cI * (TMP12)) + (TMP3 * (-cI * 
  (TMP7) + cI * (TMP8)) + TMP5 * (+cI * (TMP9) - cI * (TMP10)))); 
}


KOKKOS_FUNCTION void VVV1P0_1(Kokkos::complex<double> V2[],
    const Kokkos::complex<double> V3[], 
    const Kokkos::complex<double> COUP, 
    const double M1, const double W1, 
    Kokkos::complex<double> V1[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP10; 
  Kokkos::complex<double> TMP7; 
  Kokkos::complex<double> TMP8; 
  Kokkos::complex<double> TMP9; 
  Kokkos::complex<double> denom; 
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
  TMP9 = (V2[2] * P1[0] - V2[3] * P1[1] - V2[4] * P1[2] - V2[5] * P1[3]); 
  TMP8 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP1 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP7 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP10 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] * 
  P1[3]) - M1 * (M1 - cI * W1)); 
  V1[2] = denom * (TMP1 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[2] * (+cI * (TMP9) - cI * (TMP10)))); 
  V1[3] = denom * (TMP1 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[3] * (+cI * (TMP9) - cI * (TMP10)))); 
  V1[4] = denom * (TMP1 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[4] * (+cI * (TMP9) - cI * (TMP10)))); 
  V1[5] = denom * (TMP1 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI * 
  (TMP7) + cI * (TMP8)) + V3[5] * (+cI * (TMP9) - cI * (TMP10)))); 
}


template <typename Stride, typename Device, typename MemSpace, typename ExecSpace>
KOKKOS_FUNCTION void calculate_wavefunctions(
    Kokkos::View<int*,Stride,Device,MemSpace> cHel,
    Kokkos::View<double**,Stride,Device,MemSpace> local_mom,
    Kokkos::View<double*,ExecSpace> cIPD,
    Kokkos::View<Kokkos::complex<double>*,ExecSpace> cIPC,
    double& matrix
    )
{
  Kokkos::complex<double> amp[159];
  // Calculate wavefunctions for all processes
  Kokkos::complex<double> w[26][6];

  vxxxxx(Kokkos::subview(local_mom,0,Kokkos::ALL), 0., cHel(0), -1, w[0]); 
  vxxxxx(Kokkos::subview(local_mom,1,Kokkos::ALL), 0., cHel(1), -1, w[1]); 
  oxxxxx(Kokkos::subview(local_mom,2,Kokkos::ALL), cIPD(0), cHel(2), +1, w[2]); 
  ixxxxx(Kokkos::subview(local_mom,3,Kokkos::ALL), cIPD(0), cHel(3), -1, w[3]); 
  vxxxxx(Kokkos::subview(local_mom,4,Kokkos::ALL), 0., cHel(4), +1, w[4]); 
  vxxxxx(Kokkos::subview(local_mom,5,Kokkos::ALL), 0., cHel(5), +1, w[5]); 

  VVV1P0_1(w[0], w[1], cIPC(0), 0., 0., w[6]);

  FFV1P0_3(w[3], w[2], cIPC(1), 0., 0., w[7]);
  // Amplitude(s) for diagram number 1
  VVVV1_0(w[6], w[7], w[4], w[5], cIPC(2),
      &amp[0]);
  VVVV3_0(w[6], w[7], w[4], w[5], cIPC(2),
      &amp[1]);
  VVVV4_0(w[6], w[7], w[4], w[5], cIPC(2),
      &amp[2]);
  VVV1P0_1(w[6], w[4], cIPC(0), 0., 0.,
      w[8]);

  // Amplitude(s) for diagram number 2
  VVV1_0(w[7], w[5], w[8], cIPC(0),
      &amp[3]);
  VVV1P0_1(w[6], w[5], cIPC(0), 0., 0.,
      w[9]);
  // Amplitude(s) for diagram number 3
  VVV1_0(w[7], w[4], w[9], cIPC(0),
      &amp[4]);
  VVV1P0_1(w[4], w[5], cIPC(0), 0., 0.,
      w[10]);
  // Amplitude(s) for diagram number 4
  VVV1_0(w[6], w[7], w[10], cIPC(0),
      &amp[5]);
  FFV1_1(w[2], w[4], cIPC(1), cIPD(0),
      cIPD(1), w[11]);
  FFV1_2(w[3], w[6], cIPC(1), cIPD(0),
      cIPD(1), w[12]);
  // Amplitude(s) for diagram number 5
  FFV1_0(w[12], w[11], w[5], cIPC(1),
      &amp[6]);
  // Amplitude(s) for diagram number 6
  FFV1_0(w[3], w[11], w[9], cIPC(1),
      &amp[7]);
  FFV1_2(w[3], w[5], cIPC(1), cIPD(0),
      cIPD(1), w[13]);
  // Amplitude(s) for diagram number 7
  FFV1_0(w[13], w[11], w[6], cIPC(1),
      &amp[8]);
  FFV1_1(w[2], w[5], cIPC(1), cIPD(0),
      cIPD(1), w[14]);
  // Amplitude(s) for diagram number 8
  FFV1_0(w[12], w[14], w[4], cIPC(1),
      &amp[9]);
  // Amplitude(s) for diagram number 9
  FFV1_0(w[3], w[14], w[8], cIPC(1),
      &amp[10]);
  FFV1_2(w[3], w[4], cIPC(1), cIPD(0),
      cIPD(1), w[15]);
  // Amplitude(s) for diagram number 10
  FFV1_0(w[15], w[14], w[6], cIPC(1),
      &amp[11]);
  FFV1_1(w[2], w[6], cIPC(1), cIPD(0),
      cIPD(1), w[16]);
  // Amplitude(s) for diagram number 11
  FFV1_0(w[15], w[16], w[5], cIPC(1),
      &amp[12]);
  // Amplitude(s) for diagram number 12
  FFV1_0(w[15], w[2], w[9], cIPC(1),
      &amp[13]);
  // Amplitude(s) for diagram number 13
  FFV1_0(w[13], w[16], w[4], cIPC(1),
      &amp[14]);
  // Amplitude(s) for diagram number 14
  FFV1_0(w[13], w[2], w[8], cIPC(1),
      &amp[15]);
  // Amplitude(s) for diagram number 15
  FFV1_0(w[3], w[16], w[10], cIPC(1),
      &amp[16]);
  // Amplitude(s) for diagram number 16
  FFV1_0(w[12], w[2], w[10], cIPC(1),
      &amp[17]);
  FFV1_1(w[2], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[12]);
  FFV1_2(w[3], w[1], cIPC(1), cIPD(0),
      cIPD(1), w[16]);
  FFV1_1(w[12], w[4], cIPC(1), cIPD(0),
      cIPD(1), w[8]);
  // Amplitude(s) for diagram number 17
  FFV1_0(w[16], w[8], w[5], cIPC(1),
      &amp[18]);
  FFV1_1(w[12], w[5], cIPC(1), cIPD(0),
      cIPD(1), w[9]);
  // Amplitude(s) for diagram number 18
  FFV1_0(w[16], w[9], w[4], cIPC(1),
      &amp[19]);
  // Amplitude(s) for diagram number 19
  FFV1_0(w[16], w[12], w[10], cIPC(1),
      &amp[20]);
  VVV1P0_1(w[1], w[4], cIPC(0), 0., 0.,
      w[6]);
  FFV1P0_3(w[3], w[12], cIPC(1), 0., 0.,
      w[17]);
  // Amplitude(s) for diagram number 20
  VVV1_0(w[6], w[5], w[17], cIPC(0),
      &amp[21]);
  // Amplitude(s) for diagram number 21
  FFV1_0(w[3], w[9], w[6], cIPC(1),
      &amp[22]);
  // Amplitude(s) for diagram number 22
  FFV1_0(w[13], w[12], w[6], cIPC(1),
      &amp[23]);
  VVV1P0_1(w[1], w[5], cIPC(0), 0., 0.,
      w[18]);
  // Amplitude(s) for diagram number 23
  VVV1_0(w[18], w[4], w[17], cIPC(0),
      &amp[24]);
  // Amplitude(s) for diagram number 24
  FFV1_0(w[3], w[8], w[18], cIPC(1),
      &amp[25]);
  // Amplitude(s) for diagram number 25
  FFV1_0(w[15], w[12], w[18], cIPC(1),
      &amp[26]);
  FFV1_1(w[12], w[1], cIPC(1), cIPD(0),
      cIPD(1), w[19]);
  // Amplitude(s) for diagram number 26
  FFV1_0(w[15], w[19], w[5], cIPC(1),
      &amp[27]);
  // Amplitude(s) for diagram number 27
  FFV1_0(w[15], w[9], w[1], cIPC(1),
      &amp[28]);
  // Amplitude(s) for diagram number 28
  FFV1_0(w[13], w[19], w[4], cIPC(1),
      &amp[29]);
  // Amplitude(s) for diagram number 29
  FFV1_0(w[13], w[8], w[1], cIPC(1),
      &amp[30]);
  // Amplitude(s) for diagram number 30
  FFV1_0(w[3], w[19], w[10], cIPC(1),
      &amp[31]);
  // Amplitude(s) for diagram number 31
  VVV1_0(w[1], w[10], w[17], cIPC(0),
      &amp[32]);
  VVVV1P0_1(w[1], w[4], w[5], cIPC(2), 0.,
      0., w[17]);
  VVVV3P0_1(w[1], w[4], w[5], cIPC(2), 0.,
      0., w[19]);
  VVVV4P0_1(w[1], w[4], w[5], cIPC(2), 0.,
      0., w[8]);
  // Amplitude(s) for diagram number 32
  FFV1_0(w[3], w[12], w[17], cIPC(1),
      &amp[33]);
  FFV1_0(w[3], w[12], w[19], cIPC(1),
      &amp[34]);
  FFV1_0(w[3], w[12], w[8], cIPC(1),
      &amp[35]);
  FFV1_2(w[3], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[12]);
  FFV1_1(w[2], w[1], cIPC(1), cIPD(0),
      cIPD(1), w[9]);
  FFV1_2(w[12], w[4], cIPC(1), cIPD(0),
      cIPD(1), w[20]);
  // Amplitude(s) for diagram number 33
  FFV1_0(w[20], w[9], w[5], cIPC(1),
      &amp[36]);
  FFV1_2(w[12], w[5], cIPC(1), cIPD(0),
      cIPD(1), w[21]);
  // Amplitude(s) for diagram number 34
  FFV1_0(w[21], w[9], w[4], cIPC(1),
      &amp[37]);
  // Amplitude(s) for diagram number 35
  FFV1_0(w[12], w[9], w[10], cIPC(1),
      &amp[38]);
  FFV1P0_3(w[12], w[2], cIPC(1), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 36
  VVV1_0(w[6], w[5], w[22], cIPC(0),
      &amp[39]);
  // Amplitude(s) for diagram number 37
  FFV1_0(w[21], w[2], w[6], cIPC(1),
      &amp[40]);
  // Amplitude(s) for diagram number 38
  FFV1_0(w[12], w[14], w[6], cIPC(1),
      &amp[41]);
  // Amplitude(s) for diagram number 39
  VVV1_0(w[18], w[4], w[22], cIPC(0),
      &amp[42]);
  // Amplitude(s) for diagram number 40
  FFV1_0(w[20], w[2], w[18], cIPC(1),
      &amp[43]);
  // Amplitude(s) for diagram number 41
  FFV1_0(w[12], w[11], w[18], cIPC(1),
      &amp[44]);
  FFV1_2(w[12], w[1], cIPC(1), cIPD(0),
      cIPD(1), w[23]);
  // Amplitude(s) for diagram number 42
  FFV1_0(w[23], w[11], w[5], cIPC(1),
      &amp[45]);
  // Amplitude(s) for diagram number 43
  FFV1_0(w[21], w[11], w[1], cIPC(1),
      &amp[46]);
  // Amplitude(s) for diagram number 44
  FFV1_0(w[23], w[14], w[4], cIPC(1),
      &amp[47]);
  // Amplitude(s) for diagram number 45
  FFV1_0(w[20], w[14], w[1], cIPC(1),
      &amp[48]);
  // Amplitude(s) for diagram number 46
  FFV1_0(w[23], w[2], w[10], cIPC(1),
      &amp[49]);
  // Amplitude(s) for diagram number 47
  VVV1_0(w[1], w[10], w[22], cIPC(0),
      &amp[50]);
  // Amplitude(s) for diagram number 48
  FFV1_0(w[12], w[2], w[17], cIPC(1),
      &amp[51]);
  FFV1_0(w[12], w[2], w[19], cIPC(1),
      &amp[52]);
  FFV1_0(w[12], w[2], w[8], cIPC(1),
      &amp[53]);
  VVV1P0_1(w[0], w[4], cIPC(0), 0., 0.,
      w[12]);
  FFV1_2(w[3], w[12], cIPC(1), cIPD(0),
      cIPD(1), w[22]);
  // Amplitude(s) for diagram number 49
  FFV1_0(w[22], w[9], w[5], cIPC(1),
      &amp[54]);
  VVV1P0_1(w[12], w[5], cIPC(0), 0., 0.,
      w[23]);
  // Amplitude(s) for diagram number 50
  FFV1_0(w[3], w[9], w[23], cIPC(1),
      &amp[55]);
  // Amplitude(s) for diagram number 51
  FFV1_0(w[13], w[9], w[12], cIPC(1),
      &amp[56]);
  FFV1_1(w[2], w[12], cIPC(1), cIPD(0),
      cIPD(1), w[20]);
  // Amplitude(s) for diagram number 52
  FFV1_0(w[16], w[20], w[5], cIPC(1),
      &amp[57]);
  // Amplitude(s) for diagram number 53
  FFV1_0(w[16], w[2], w[23], cIPC(1),
      &amp[58]);
  // Amplitude(s) for diagram number 54
  FFV1_0(w[16], w[14], w[12], cIPC(1),
      &amp[59]);
  // Amplitude(s) for diagram number 55
  FFV1_0(w[3], w[20], w[18], cIPC(1),
      &amp[60]);
  // Amplitude(s) for diagram number 56
  FFV1_0(w[22], w[2], w[18], cIPC(1),
      &amp[61]);
  // Amplitude(s) for diagram number 57
  VVV1_0(w[12], w[18], w[7], cIPC(0),
      &amp[62]);
  // Amplitude(s) for diagram number 58
  VVVV1_0(w[12], w[1], w[7], w[5], cIPC(2),
      &amp[63]);
  VVVV3_0(w[12], w[1], w[7], w[5], cIPC(2),
      &amp[64]);
  VVVV4_0(w[12], w[1], w[7], w[5], cIPC(2),
      &amp[65]);
  VVV1P0_1(w[12], w[1], cIPC(0), 0., 0.,
      w[21]);
  // Amplitude(s) for diagram number 59
  VVV1_0(w[7], w[5], w[21], cIPC(0),
      &amp[66]);
  // Amplitude(s) for diagram number 60
  VVV1_0(w[1], w[7], w[23], cIPC(0),
      &amp[67]);
  // Amplitude(s) for diagram number 61
  FFV1_0(w[3], w[14], w[21], cIPC(1),
      &amp[68]);
  // Amplitude(s) for diagram number 62
  FFV1_0(w[22], w[14], w[1], cIPC(1),
      &amp[69]);
  // Amplitude(s) for diagram number 63
  FFV1_0(w[13], w[2], w[21], cIPC(1),
      &amp[70]);
  // Amplitude(s) for diagram number 64
  FFV1_0(w[13], w[20], w[1], cIPC(1),
      &amp[71]);
  VVV1P0_1(w[0], w[5], cIPC(0), 0., 0.,
      w[20]);
  FFV1_2(w[3], w[20], cIPC(1), cIPD(0),
      cIPD(1), w[21]);
  // Amplitude(s) for diagram number 65
  FFV1_0(w[21], w[9], w[4], cIPC(1),
      &amp[72]);
  VVV1P0_1(w[20], w[4], cIPC(0), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 66
  FFV1_0(w[3], w[9], w[22], cIPC(1),
      &amp[73]);
  // Amplitude(s) for diagram number 67
  FFV1_0(w[15], w[9], w[20], cIPC(1),
      &amp[74]);
  FFV1_1(w[2], w[20], cIPC(1), cIPD(0),
      cIPD(1), w[23]);
  // Amplitude(s) for diagram number 68
  FFV1_0(w[16], w[23], w[4], cIPC(1),
      &amp[75]);
  // Amplitude(s) for diagram number 69
  FFV1_0(w[16], w[2], w[22], cIPC(1),
      &amp[76]);
  // Amplitude(s) for diagram number 70
  FFV1_0(w[16], w[11], w[20], cIPC(1),
      &amp[77]);
  // Amplitude(s) for diagram number 71
  FFV1_0(w[3], w[23], w[6], cIPC(1),
      &amp[78]);
  // Amplitude(s) for diagram number 72
  FFV1_0(w[21], w[2], w[6], cIPC(1),
      &amp[79]);
  // Amplitude(s) for diagram number 73
  VVV1_0(w[20], w[6], w[7], cIPC(0),
      &amp[80]);
  // Amplitude(s) for diagram number 74
  VVVV1_0(w[20], w[1], w[7], w[4], cIPC(2),
      &amp[81]);
  VVVV3_0(w[20], w[1], w[7], w[4], cIPC(2),
      &amp[82]);
  VVVV4_0(w[20], w[1], w[7], w[4], cIPC(2),
      &amp[83]);
  VVV1P0_1(w[20], w[1], cIPC(0), 0., 0.,
      w[12]);
  // Amplitude(s) for diagram number 75
  VVV1_0(w[7], w[4], w[12], cIPC(0),
      &amp[84]);
  // Amplitude(s) for diagram number 76
  VVV1_0(w[1], w[7], w[22], cIPC(0),
      &amp[85]);
  // Amplitude(s) for diagram number 77
  FFV1_0(w[3], w[11], w[12], cIPC(1),
      &amp[86]);
  // Amplitude(s) for diagram number 78
  FFV1_0(w[21], w[11], w[1], cIPC(1),
      &amp[87]);
  // Amplitude(s) for diagram number 79
  FFV1_0(w[15], w[2], w[12], cIPC(1),
      &amp[88]);
  // Amplitude(s) for diagram number 80
  FFV1_0(w[15], w[23], w[1], cIPC(1),
      &amp[89]);
  FFV1_1(w[9], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[23]);
  // Amplitude(s) for diagram number 81
  FFV1_0(w[15], w[23], w[5], cIPC(1),
      &amp[90]);
  FFV1_2(w[15], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[12]);
  // Amplitude(s) for diagram number 82
  FFV1_0(w[12], w[9], w[5], cIPC(1),
      &amp[91]);
  // Amplitude(s) for diagram number 83
  FFV1_0(w[13], w[23], w[4], cIPC(1),
      &amp[92]);
  FFV1_2(w[13], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[21]);
  // Amplitude(s) for diagram number 84
  FFV1_0(w[21], w[9], w[4], cIPC(1),
      &amp[93]);
  // Amplitude(s) for diagram number 85
  FFV1_0(w[3], w[23], w[10], cIPC(1),
      &amp[94]);
  VVV1P0_1(w[0], w[10], cIPC(0), 0., 0.,
      w[23]);
  // Amplitude(s) for diagram number 86
  FFV1_0(w[3], w[9], w[23], cIPC(1),
      &amp[95]);
  FFV1_2(w[16], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[22]);
  // Amplitude(s) for diagram number 87
  FFV1_0(w[22], w[11], w[5], cIPC(1),
      &amp[96]);
  FFV1_1(w[11], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[20]);
  // Amplitude(s) for diagram number 88
  FFV1_0(w[16], w[20], w[5], cIPC(1),
      &amp[97]);
  // Amplitude(s) for diagram number 89
  FFV1_0(w[22], w[14], w[4], cIPC(1),
      &amp[98]);
  FFV1_1(w[14], w[0], cIPC(1), cIPD(0),
      cIPD(1), w[24]);
  // Amplitude(s) for diagram number 90
  FFV1_0(w[16], w[24], w[4], cIPC(1),
      &amp[99]);
  // Amplitude(s) for diagram number 91
  FFV1_0(w[22], w[2], w[10], cIPC(1),
      &amp[100]);
  // Amplitude(s) for diagram number 92
  FFV1_0(w[16], w[2], w[23], cIPC(1),
      &amp[101]);
  // Amplitude(s) for diagram number 93
  VVVV1_0(w[0], w[6], w[7], w[5], cIPC(2),
      &amp[102]);
  VVVV3_0(w[0], w[6], w[7], w[5], cIPC(2),
      &amp[103]);
  VVVV4_0(w[0], w[6], w[7], w[5], cIPC(2),
      &amp[104]);
  VVV1P0_1(w[0], w[6], cIPC(0), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 94
  VVV1_0(w[7], w[5], w[22], cIPC(0),
      &amp[105]);
  VVV1P0_1(w[0], w[7], cIPC(0), 0., 0.,
      w[25]);
  // Amplitude(s) for diagram number 95
  VVV1_0(w[6], w[5], w[25], cIPC(0),
      &amp[106]);
  // Amplitude(s) for diagram number 96
  FFV1_0(w[3], w[14], w[22], cIPC(1),
      &amp[107]);
  // Amplitude(s) for diagram number 97
  FFV1_0(w[3], w[24], w[6], cIPC(1),
      &amp[108]);
  // Amplitude(s) for diagram number 98
  FFV1_0(w[13], w[2], w[22], cIPC(1),
      &amp[109]);
  // Amplitude(s) for diagram number 99
  FFV1_0(w[21], w[2], w[6], cIPC(1),
      &amp[110]);
  // Amplitude(s) for diagram number 100
  VVVV1_0(w[0], w[18], w[7], w[4], cIPC(2),
      &amp[111]);
  VVVV3_0(w[0], w[18], w[7], w[4], cIPC(2),
      &amp[112]);
  VVVV4_0(w[0], w[18], w[7], w[4], cIPC(2),
      &amp[113]);
  VVV1P0_1(w[0], w[18], cIPC(0), 0., 0.,
      w[6]);
  // Amplitude(s) for diagram number 101
  VVV1_0(w[7], w[4], w[6], cIPC(0),
      &amp[114]);
  // Amplitude(s) for diagram number 102
  VVV1_0(w[18], w[4], w[25], cIPC(0),
      &amp[115]);
  // Amplitude(s) for diagram number 103
  FFV1_0(w[3], w[11], w[6], cIPC(1),
      &amp[116]);
  // Amplitude(s) for diagram number 104
  FFV1_0(w[3], w[20], w[18], cIPC(1),
      &amp[117]);
  // Amplitude(s) for diagram number 105
  FFV1_0(w[15], w[2], w[6], cIPC(1),
      &amp[118]);
  // Amplitude(s) for diagram number 106
  FFV1_0(w[12], w[2], w[18], cIPC(1),
      &amp[119]);
  // Amplitude(s) for diagram number 107
  VVVV1_0(w[0], w[1], w[7], w[10], cIPC(2),
      &amp[120]);
  VVVV3_0(w[0], w[1], w[7], w[10], cIPC(2),
      &amp[121]);
  VVVV4_0(w[0], w[1], w[7], w[10], cIPC(2),
      &amp[122]);
  // Amplitude(s) for diagram number 108
  VVV1_0(w[1], w[10], w[25], cIPC(0),
      &amp[123]);
  // Amplitude(s) for diagram number 109
  VVV1_0(w[1], w[7], w[23], cIPC(0),
      &amp[124]);
  // Amplitude(s) for diagram number 110
  FFV1_0(w[13], w[20], w[1], cIPC(1),
      &amp[125]);
  // Amplitude(s) for diagram number 111
  FFV1_0(w[21], w[11], w[1], cIPC(1),
      &amp[126]);
  // Amplitude(s) for diagram number 112
  FFV1_0(w[15], w[24], w[1], cIPC(1),
      &amp[127]);
  // Amplitude(s) for diagram number 113
  FFV1_0(w[12], w[14], w[1], cIPC(1),
      &amp[128]);
  VVVV1P0_1(w[0], w[1], w[4], cIPC(2), 0.,
      0., w[12]);
  VVVV3P0_1(w[0], w[1], w[4], cIPC(2), 0.,
      0., w[24]);
  VVVV4P0_1(w[0], w[1], w[4], cIPC(2), 0.,
      0., w[21]);
  // Amplitude(s) for diagram number 114
  VVV1_0(w[12], w[7], w[5], cIPC(0),
      &amp[129]);
  VVV1_0(w[24], w[7], w[5], cIPC(0),
      &amp[130]);
  VVV1_0(w[21], w[7], w[5], cIPC(0),
      &amp[131]);
  // Amplitude(s) for diagram number 115
  FFV1_0(w[3], w[14], w[12], cIPC(1),
      &amp[132]);
  FFV1_0(w[3], w[14], w[24], cIPC(1),
      &amp[133]);
  FFV1_0(w[3], w[14], w[21], cIPC(1),
      &amp[134]);
  // Amplitude(s) for diagram number 116
  FFV1_0(w[13], w[2], w[12], cIPC(1),
      &amp[135]);
  FFV1_0(w[13], w[2], w[24], cIPC(1),
      &amp[136]);
  FFV1_0(w[13], w[2], w[21], cIPC(1),
      &amp[137]);
  VVVV1P0_1(w[0], w[1], w[5], cIPC(2), 0.,
      0., w[21]);
  VVVV3P0_1(w[0], w[1], w[5], cIPC(2), 0.,
      0., w[13]);
  VVVV4P0_1(w[0], w[1], w[5], cIPC(2), 0.,
      0., w[24]);
  // Amplitude(s) for diagram number 117
  VVV1_0(w[21], w[7], w[4], cIPC(0),
      &amp[138]);
  VVV1_0(w[13], w[7], w[4], cIPC(0),
      &amp[139]);
  VVV1_0(w[24], w[7], w[4], cIPC(0),
      &amp[140]);
  // Amplitude(s) for diagram number 118
  FFV1_0(w[3], w[11], w[21], cIPC(1),
      &amp[141]);
  FFV1_0(w[3], w[11], w[13], cIPC(1),
      &amp[142]);
  FFV1_0(w[3], w[11], w[24], cIPC(1),
      &amp[143]);
  // Amplitude(s) for diagram number 119
  FFV1_0(w[15], w[2], w[21], cIPC(1),
      &amp[144]);
  FFV1_0(w[15], w[2], w[13], cIPC(1),
      &amp[145]);
  FFV1_0(w[15], w[2], w[24], cIPC(1),
      &amp[146]);
  VVVV1P0_1(w[0], w[4], w[5], cIPC(2), 0.,
      0., w[24]);
  VVVV3P0_1(w[0], w[4], w[5], cIPC(2), 0.,
      0., w[15]);
  VVVV4P0_1(w[0], w[4], w[5], cIPC(2), 0.,
      0., w[13]);
  // Amplitude(s) for diagram number 120
  FFV1_0(w[3], w[9], w[24], cIPC(1),
      &amp[147]);
  FFV1_0(w[3], w[9], w[15], cIPC(1),
      &amp[148]);
  FFV1_0(w[3], w[9], w[13], cIPC(1),
      &amp[149]);
  // Amplitude(s) for diagram number 121
  FFV1_0(w[16], w[2], w[24], cIPC(1),
      &amp[150]);
  FFV1_0(w[16], w[2], w[15], cIPC(1),
      &amp[151]);
  FFV1_0(w[16], w[2], w[13], cIPC(1),
      &amp[152]);
  // Amplitude(s) for diagram number 122
  VVV1_0(w[24], w[1], w[7], cIPC(0),
      &amp[153]);
  VVV1_0(w[15], w[1], w[7], cIPC(0),
      &amp[154]);
  VVV1_0(w[13], w[1], w[7], cIPC(0),
      &amp[155]);
  // Amplitude(s) for diagram number 123
  VVV1_0(w[0], w[17], w[7], cIPC(0),
      &amp[156]);
  VVV1_0(w[0], w[19], w[7], cIPC(0),
      &amp[157]);
  VVV1_0(w[0], w[8], w[7], cIPC(0),
      &amp[158]);

  // double CPPProcess::matrix_1_gg_ttxgg() {
  int i, j; 
  // Local variables

  // const int ngraphs = 2;
  const int ncolor = 24; 
  Kokkos::complex<double> ztemp; 
  Kokkos::complex<double> jamp[ncolor]; 
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

  // Calculate color flows
  jamp[0] = +Kokkos::complex<double> (0, 1) * amp[0] + Kokkos::complex<double>
      (0, 1) * amp[1] + Kokkos::complex<double> (0, 1) * amp[3] +
      Kokkos::complex<double> (0, 1) * amp[5] + Kokkos::complex<double> (0, 1)
      * amp[14] + amp[15] + amp[16] + amp[21] + Kokkos::complex<double> (0, 1)
      * amp[23] - amp[29] + Kokkos::complex<double> (0, 1) * amp[31] + amp[32]
      + amp[33] - amp[35] + Kokkos::complex<double> (0, 1) * amp[102] +
      Kokkos::complex<double> (0, 1) * amp[103] + Kokkos::complex<double> (0,
      1) * amp[105] + Kokkos::complex<double> (0, 1) * amp[106] + amp[109] +
      Kokkos::complex<double> (0, 1) * amp[120] + Kokkos::complex<double> (0,
      1) * amp[121] + Kokkos::complex<double> (0, 1) * amp[123] +
      Kokkos::complex<double> (0, 1) * amp[129] - Kokkos::complex<double> (0,
      1) * amp[131] + amp[135] - amp[137] - Kokkos::complex<double> (0, 1) *
      amp[156] + Kokkos::complex<double> (0, 1) * amp[158];
  jamp[1] = -Kokkos::complex<double> (0, 1) * amp[0] + Kokkos::complex<double>
      (0, 1) * amp[2] + Kokkos::complex<double> (0, 1) * amp[4] -
      Kokkos::complex<double> (0, 1) * amp[5] + Kokkos::complex<double> (0, 1)
      * amp[12] + amp[13] - amp[16] + amp[24] + Kokkos::complex<double> (0, 1)
      * amp[26] - amp[27] - Kokkos::complex<double> (0, 1) * amp[31] - amp[32]
      - amp[33] - amp[34] + Kokkos::complex<double> (0, 1) * amp[111] +
      Kokkos::complex<double> (0, 1) * amp[112] + Kokkos::complex<double> (0,
      1) * amp[114] + Kokkos::complex<double> (0, 1) * amp[115] + amp[118] -
      Kokkos::complex<double> (0, 1) * amp[120] - Kokkos::complex<double> (0,
      1) * amp[121] - Kokkos::complex<double> (0, 1) * amp[123] +
      Kokkos::complex<double> (0, 1) * amp[138] - Kokkos::complex<double> (0,
      1) * amp[140] + amp[144] - amp[146] + Kokkos::complex<double> (0, 1) *
      amp[156] + Kokkos::complex<double> (0, 1) * amp[157];
  jamp[2] = -amp[21] - Kokkos::complex<double> (0, 1) * amp[23] - amp[24] +
      Kokkos::complex<double> (0, 1) * amp[25] - amp[30] + amp[34] + amp[35] +
      amp[60] - Kokkos::complex<double> (0, 1) * amp[62] +
      Kokkos::complex<double> (0, 1) * amp[63] + Kokkos::complex<double> (0, 1)
      * amp[64] + Kokkos::complex<double> (0, 1) * amp[66] + amp[70] +
      Kokkos::complex<double> (0, 1) * amp[71] - Kokkos::complex<double> (0, 1)
      * amp[102] - Kokkos::complex<double> (0, 1) * amp[103] -
      Kokkos::complex<double> (0, 1) * amp[105] - Kokkos::complex<double> (0,
      1) * amp[106] - amp[109] - Kokkos::complex<double> (0, 1) * amp[112] -
      Kokkos::complex<double> (0, 1) * amp[113] - Kokkos::complex<double> (0,
      1) * amp[115] - Kokkos::complex<double> (0, 1) * amp[129] -
      Kokkos::complex<double> (0, 1) * amp[130] - amp[135] - amp[136] -
      Kokkos::complex<double> (0, 1) * amp[157] - Kokkos::complex<double> (0,
      1) * amp[158];
  jamp[3] = -amp[18] + Kokkos::complex<double> (0, 1) * amp[20] + amp[24] -
      Kokkos::complex<double> (0, 1) * amp[25] - amp[32] - amp[33] - amp[34] +
      Kokkos::complex<double> (0, 1) * amp[57] + amp[58] - amp[60] +
      Kokkos::complex<double> (0, 1) * amp[62] - Kokkos::complex<double> (0, 1)
      * amp[64] - Kokkos::complex<double> (0, 1) * amp[65] -
      Kokkos::complex<double> (0, 1) * amp[67] + amp[101] +
      Kokkos::complex<double> (0, 1) * amp[112] + Kokkos::complex<double> (0,
      1) * amp[113] + Kokkos::complex<double> (0, 1) * amp[115] -
      Kokkos::complex<double> (0, 1) * amp[121] - Kokkos::complex<double> (0,
      1) * amp[122] - Kokkos::complex<double> (0, 1) * amp[123] -
      Kokkos::complex<double> (0, 1) * amp[124] + amp[150] - amp[152] -
      Kokkos::complex<double> (0, 1) * amp[153] + Kokkos::complex<double> (0,
      1) * amp[155] + Kokkos::complex<double> (0, 1) * amp[156] +
      Kokkos::complex<double> (0, 1) * amp[157];
  jamp[4] = -amp[21] + Kokkos::complex<double> (0, 1) * amp[22] - amp[24] -
      Kokkos::complex<double> (0, 1) * amp[26] - amp[28] + amp[34] + amp[35] +
      amp[78] - Kokkos::complex<double> (0, 1) * amp[80] +
      Kokkos::complex<double> (0, 1) * amp[81] + Kokkos::complex<double> (0, 1)
      * amp[82] + Kokkos::complex<double> (0, 1) * amp[84] + amp[88] +
      Kokkos::complex<double> (0, 1) * amp[89] - Kokkos::complex<double> (0, 1)
      * amp[103] - Kokkos::complex<double> (0, 1) * amp[104] -
      Kokkos::complex<double> (0, 1) * amp[106] - Kokkos::complex<double> (0,
      1) * amp[111] - Kokkos::complex<double> (0, 1) * amp[112] -
      Kokkos::complex<double> (0, 1) * amp[114] - Kokkos::complex<double> (0,
      1) * amp[115] - amp[118] - Kokkos::complex<double> (0, 1) * amp[138] -
      Kokkos::complex<double> (0, 1) * amp[139] - amp[144] - amp[145] -
      Kokkos::complex<double> (0, 1) * amp[157] - Kokkos::complex<double> (0,
      1) * amp[158];
  jamp[5] = -amp[19] - Kokkos::complex<double> (0, 1) * amp[20] + amp[21] -
      Kokkos::complex<double> (0, 1) * amp[22] + amp[32] + amp[33] - amp[35] +
      Kokkos::complex<double> (0, 1) * amp[75] + amp[76] - amp[78] +
      Kokkos::complex<double> (0, 1) * amp[80] - Kokkos::complex<double> (0, 1)
      * amp[82] - Kokkos::complex<double> (0, 1) * amp[83] -
      Kokkos::complex<double> (0, 1) * amp[85] - amp[101] +
      Kokkos::complex<double> (0, 1) * amp[103] + Kokkos::complex<double> (0,
      1) * amp[104] + Kokkos::complex<double> (0, 1) * amp[106] +
      Kokkos::complex<double> (0, 1) * amp[121] + Kokkos::complex<double> (0,
      1) * amp[122] + Kokkos::complex<double> (0, 1) * amp[123] +
      Kokkos::complex<double> (0, 1) * amp[124] - amp[150] - amp[151] +
      Kokkos::complex<double> (0, 1) * amp[153] + Kokkos::complex<double> (0,
      1) * amp[154] - Kokkos::complex<double> (0, 1) * amp[156] +
      Kokkos::complex<double> (0, 1) * amp[158];
  jamp[6] = -Kokkos::complex<double> (0, 1) * amp[0] - Kokkos::complex<double>
      (0, 1) * amp[1] - Kokkos::complex<double> (0, 1) * amp[3] -
      Kokkos::complex<double> (0, 1) * amp[5] - Kokkos::complex<double> (0, 1)
      * amp[14] - amp[15] - amp[16] + amp[55] + Kokkos::complex<double> (0, 1)
      * amp[56] - Kokkos::complex<double> (0, 1) * amp[63] +
      Kokkos::complex<double> (0, 1) * amp[65] - Kokkos::complex<double> (0, 1)
      * amp[66] + Kokkos::complex<double> (0, 1) * amp[67] - amp[70] - amp[92]
      + Kokkos::complex<double> (0, 1) * amp[94] + amp[95] -
      Kokkos::complex<double> (0, 1) * amp[120] + Kokkos::complex<double> (0,
      1) * amp[122] + Kokkos::complex<double> (0, 1) * amp[124] +
      Kokkos::complex<double> (0, 1) * amp[130] + Kokkos::complex<double> (0,
      1) * amp[131] + amp[136] + amp[137] + amp[147] - amp[149] +
      Kokkos::complex<double> (0, 1) * amp[153] - Kokkos::complex<double> (0,
      1) * amp[155];
  jamp[7] = +Kokkos::complex<double> (0, 1) * amp[0] - Kokkos::complex<double>
      (0, 1) * amp[2] - Kokkos::complex<double> (0, 1) * amp[4] +
      Kokkos::complex<double> (0, 1) * amp[5] - Kokkos::complex<double> (0, 1)
      * amp[12] - amp[13] + amp[16] + amp[73] + Kokkos::complex<double> (0, 1)
      * amp[74] - Kokkos::complex<double> (0, 1) * amp[81] +
      Kokkos::complex<double> (0, 1) * amp[83] - Kokkos::complex<double> (0, 1)
      * amp[84] + Kokkos::complex<double> (0, 1) * amp[85] - amp[88] - amp[90]
      - Kokkos::complex<double> (0, 1) * amp[94] - amp[95] +
      Kokkos::complex<double> (0, 1) * amp[120] - Kokkos::complex<double> (0,
      1) * amp[122] - Kokkos::complex<double> (0, 1) * amp[124] +
      Kokkos::complex<double> (0, 1) * amp[139] + Kokkos::complex<double> (0,
      1) * amp[140] + amp[145] + amp[146] - amp[147] - amp[148] -
      Kokkos::complex<double> (0, 1) * amp[153] - Kokkos::complex<double> (0,
      1) * amp[154];
  jamp[8] = -amp[55] - Kokkos::complex<double> (0, 1) * amp[56] +
      Kokkos::complex<double> (0, 1) * amp[63] - Kokkos::complex<double> (0, 1)
      * amp[65] + Kokkos::complex<double> (0, 1) * amp[66] -
      Kokkos::complex<double> (0, 1) * amp[67] + amp[70] +
      Kokkos::complex<double> (0, 1) * amp[72] - amp[73] + amp[79] +
      Kokkos::complex<double> (0, 1) * amp[80] - Kokkos::complex<double> (0, 1)
      * amp[82] - Kokkos::complex<double> (0, 1) * amp[83] -
      Kokkos::complex<double> (0, 1) * amp[85] - amp[93] -
      Kokkos::complex<double> (0, 1) * amp[102] + Kokkos::complex<double> (0,
      1) * amp[104] - Kokkos::complex<double> (0, 1) * amp[105] - amp[109] +
      Kokkos::complex<double> (0, 1) * amp[110] - Kokkos::complex<double> (0,
      1) * amp[129] - Kokkos::complex<double> (0, 1) * amp[130] - amp[135] -
      amp[136] + amp[148] + amp[149] + Kokkos::complex<double> (0, 1) *
      amp[154] + Kokkos::complex<double> (0, 1) * amp[155];
  jamp[9] = -amp[37] + Kokkos::complex<double> (0, 1) * amp[38] + amp[39] +
      Kokkos::complex<double> (0, 1) * amp[40] + amp[50] + amp[51] - amp[53] -
      Kokkos::complex<double> (0, 1) * amp[72] + amp[73] - amp[79] -
      Kokkos::complex<double> (0, 1) * amp[80] + Kokkos::complex<double> (0, 1)
      * amp[82] + Kokkos::complex<double> (0, 1) * amp[83] +
      Kokkos::complex<double> (0, 1) * amp[85] - amp[95] -
      Kokkos::complex<double> (0, 1) * amp[103] - Kokkos::complex<double> (0,
      1) * amp[104] - Kokkos::complex<double> (0, 1) * amp[106] -
      Kokkos::complex<double> (0, 1) * amp[121] - Kokkos::complex<double> (0,
      1) * amp[122] - Kokkos::complex<double> (0, 1) * amp[123] -
      Kokkos::complex<double> (0, 1) * amp[124] - amp[147] - amp[148] -
      Kokkos::complex<double> (0, 1) * amp[153] - Kokkos::complex<double> (0,
      1) * amp[154] + Kokkos::complex<double> (0, 1) * amp[156] -
      Kokkos::complex<double> (0, 1) * amp[158];
  jamp[10] = +Kokkos::complex<double> (0, 1) * amp[54] - amp[55] + amp[61] +
      Kokkos::complex<double> (0, 1) * amp[62] - Kokkos::complex<double> (0, 1)
      * amp[64] - Kokkos::complex<double> (0, 1) * amp[65] -
      Kokkos::complex<double> (0, 1) * amp[67] - amp[73] -
      Kokkos::complex<double> (0, 1) * amp[74] + Kokkos::complex<double> (0, 1)
      * amp[81] - Kokkos::complex<double> (0, 1) * amp[83] +
      Kokkos::complex<double> (0, 1) * amp[84] - Kokkos::complex<double> (0, 1)
      * amp[85] + amp[88] - amp[91] - Kokkos::complex<double> (0, 1) * amp[111]
      + Kokkos::complex<double> (0, 1) * amp[113] - Kokkos::complex<double> (0,
      1) * amp[114] - amp[118] + Kokkos::complex<double> (0, 1) * amp[119] -
      Kokkos::complex<double> (0, 1) * amp[138] - Kokkos::complex<double> (0,
      1) * amp[139] - amp[144] - amp[145] + amp[148] + amp[149] +
      Kokkos::complex<double> (0, 1) * amp[154] + Kokkos::complex<double> (0,
      1) * amp[155];
  jamp[11] = -amp[36] - Kokkos::complex<double> (0, 1) * amp[38] + amp[42] +
      Kokkos::complex<double> (0, 1) * amp[43] - amp[50] - amp[51] - amp[52] -
      Kokkos::complex<double> (0, 1) * amp[54] + amp[55] - amp[61] -
      Kokkos::complex<double> (0, 1) * amp[62] + Kokkos::complex<double> (0, 1)
      * amp[64] + Kokkos::complex<double> (0, 1) * amp[65] +
      Kokkos::complex<double> (0, 1) * amp[67] + amp[95] -
      Kokkos::complex<double> (0, 1) * amp[112] - Kokkos::complex<double> (0,
      1) * amp[113] - Kokkos::complex<double> (0, 1) * amp[115] +
      Kokkos::complex<double> (0, 1) * amp[121] + Kokkos::complex<double> (0,
      1) * amp[122] + Kokkos::complex<double> (0, 1) * amp[123] +
      Kokkos::complex<double> (0, 1) * amp[124] + amp[147] - amp[149] +
      Kokkos::complex<double> (0, 1) * amp[153] - Kokkos::complex<double> (0,
      1) * amp[155] - Kokkos::complex<double> (0, 1) * amp[156] -
      Kokkos::complex<double> (0, 1) * amp[157];
  jamp[12] = -Kokkos::complex<double> (0, 1) * amp[1] - Kokkos::complex<double>
      (0, 1) * amp[2] - Kokkos::complex<double> (0, 1) * amp[3] -
      Kokkos::complex<double> (0, 1) * amp[4] + amp[7] +
      Kokkos::complex<double> (0, 1) * amp[8] - amp[15] - amp[60] +
      Kokkos::complex<double> (0, 1) * amp[62] - Kokkos::complex<double> (0, 1)
      * amp[63] - Kokkos::complex<double> (0, 1) * amp[64] -
      Kokkos::complex<double> (0, 1) * amp[66] - amp[70] -
      Kokkos::complex<double> (0, 1) * amp[71] - Kokkos::complex<double> (0, 1)
      * amp[111] + Kokkos::complex<double> (0, 1) * amp[113] -
      Kokkos::complex<double> (0, 1) * amp[114] + amp[116] +
      Kokkos::complex<double> (0, 1) * amp[117] - amp[125] +
      Kokkos::complex<double> (0, 1) * amp[130] + Kokkos::complex<double> (0,
      1) * amp[131] + amp[136] + amp[137] - Kokkos::complex<double> (0, 1) *
      amp[138] + Kokkos::complex<double> (0, 1) * amp[140] + amp[141] -
      amp[143];
  jamp[13] = -Kokkos::complex<double> (0, 1) * amp[57] - amp[58] + amp[60] -
      Kokkos::complex<double> (0, 1) * amp[62] + Kokkos::complex<double> (0, 1)
      * amp[64] + Kokkos::complex<double> (0, 1) * amp[65] +
      Kokkos::complex<double> (0, 1) * amp[67] - amp[76] +
      Kokkos::complex<double> (0, 1) * amp[77] - Kokkos::complex<double> (0, 1)
      * amp[81] + Kokkos::complex<double> (0, 1) * amp[83] -
      Kokkos::complex<double> (0, 1) * amp[84] + Kokkos::complex<double> (0, 1)
      * amp[85] + amp[86] - amp[97] + Kokkos::complex<double> (0, 1) * amp[111]
      - Kokkos::complex<double> (0, 1) * amp[113] + Kokkos::complex<double> (0,
      1) * amp[114] - amp[116] - Kokkos::complex<double> (0, 1) * amp[117] +
      Kokkos::complex<double> (0, 1) * amp[138] + Kokkos::complex<double> (0,
      1) * amp[139] - amp[141] - amp[142] + amp[151] + amp[152] -
      Kokkos::complex<double> (0, 1) * amp[154] - Kokkos::complex<double> (0,
      1) * amp[155];
  jamp[14] = +Kokkos::complex<double> (0, 1) * amp[1] + Kokkos::complex<double>
      (0, 1) * amp[2] + Kokkos::complex<double> (0, 1) * amp[3] +
      Kokkos::complex<double> (0, 1) * amp[4] - amp[7] -
      Kokkos::complex<double> (0, 1) * amp[8] + amp[15] - amp[79] -
      Kokkos::complex<double> (0, 1) * amp[80] + Kokkos::complex<double> (0, 1)
      * amp[81] + Kokkos::complex<double> (0, 1) * amp[82] +
      Kokkos::complex<double> (0, 1) * amp[84] - amp[86] +
      Kokkos::complex<double> (0, 1) * amp[87] + Kokkos::complex<double> (0, 1)
      * amp[102] - Kokkos::complex<double> (0, 1) * amp[104] +
      Kokkos::complex<double> (0, 1) * amp[105] + amp[109] -
      Kokkos::complex<double> (0, 1) * amp[110] - amp[126] +
      Kokkos::complex<double> (0, 1) * amp[129] - Kokkos::complex<double> (0,
      1) * amp[131] + amp[135] - amp[137] - Kokkos::complex<double> (0, 1) *
      amp[139] - Kokkos::complex<double> (0, 1) * amp[140] + amp[142] +
      amp[143];
  jamp[15] = -amp[39] - Kokkos::complex<double> (0, 1) * amp[40] - amp[42] +
      Kokkos::complex<double> (0, 1) * amp[44] - amp[46] + amp[52] + amp[53] +
      amp[79] + Kokkos::complex<double> (0, 1) * amp[80] -
      Kokkos::complex<double> (0, 1) * amp[81] - Kokkos::complex<double> (0, 1)
      * amp[82] - Kokkos::complex<double> (0, 1) * amp[84] + amp[86] -
      Kokkos::complex<double> (0, 1) * amp[87] + Kokkos::complex<double> (0, 1)
      * amp[103] + Kokkos::complex<double> (0, 1) * amp[104] +
      Kokkos::complex<double> (0, 1) * amp[106] + Kokkos::complex<double> (0,
      1) * amp[111] + Kokkos::complex<double> (0, 1) * amp[112] +
      Kokkos::complex<double> (0, 1) * amp[114] + Kokkos::complex<double> (0,
      1) * amp[115] - amp[116] + Kokkos::complex<double> (0, 1) * amp[138] +
      Kokkos::complex<double> (0, 1) * amp[139] - amp[141] - amp[142] +
      Kokkos::complex<double> (0, 1) * amp[157] + Kokkos::complex<double> (0,
      1) * amp[158];
  jamp[16] = -Kokkos::complex<double> (0, 1) * amp[0] + Kokkos::complex<double>
      (0, 1) * amp[2] + Kokkos::complex<double> (0, 1) * amp[4] -
      Kokkos::complex<double> (0, 1) * amp[5] + Kokkos::complex<double> (0, 1)
      * amp[6] - amp[7] + amp[17] + amp[76] - Kokkos::complex<double> (0, 1) *
      amp[77] + Kokkos::complex<double> (0, 1) * amp[81] -
      Kokkos::complex<double> (0, 1) * amp[83] + Kokkos::complex<double> (0, 1)
      * amp[84] - Kokkos::complex<double> (0, 1) * amp[85] - amp[86] - amp[96]
      + Kokkos::complex<double> (0, 1) * amp[100] - amp[101] -
      Kokkos::complex<double> (0, 1) * amp[120] + Kokkos::complex<double> (0,
      1) * amp[122] + Kokkos::complex<double> (0, 1) * amp[124] -
      Kokkos::complex<double> (0, 1) * amp[139] - Kokkos::complex<double> (0,
      1) * amp[140] + amp[142] + amp[143] - amp[150] - amp[151] +
      Kokkos::complex<double> (0, 1) * amp[153] + Kokkos::complex<double> (0,
      1) * amp[154];
  jamp[17] = +Kokkos::complex<double> (0, 1) * amp[0] - Kokkos::complex<double>
      (0, 1) * amp[2] - Kokkos::complex<double> (0, 1) * amp[4] +
      Kokkos::complex<double> (0, 1) * amp[5] - Kokkos::complex<double> (0, 1)
      * amp[6] + amp[7] - amp[17] + amp[42] - Kokkos::complex<double> (0, 1) *
      amp[44] - amp[45] + Kokkos::complex<double> (0, 1) * amp[49] - amp[50] -
      amp[51] - amp[52] - Kokkos::complex<double> (0, 1) * amp[111] -
      Kokkos::complex<double> (0, 1) * amp[112] - Kokkos::complex<double> (0,
      1) * amp[114] - Kokkos::complex<double> (0, 1) * amp[115] + amp[116] +
      Kokkos::complex<double> (0, 1) * amp[120] + Kokkos::complex<double> (0,
      1) * amp[121] + Kokkos::complex<double> (0, 1) * amp[123] -
      Kokkos::complex<double> (0, 1) * amp[138] + Kokkos::complex<double> (0,
      1) * amp[140] + amp[141] - amp[143] - Kokkos::complex<double> (0, 1) *
      amp[156] - Kokkos::complex<double> (0, 1) * amp[157];
  jamp[18] = -Kokkos::complex<double> (0, 1) * amp[1] - Kokkos::complex<double>
      (0, 1) * amp[2] - Kokkos::complex<double> (0, 1) * amp[3] -
      Kokkos::complex<double> (0, 1) * amp[4] + amp[10] +
      Kokkos::complex<double> (0, 1) * amp[11] - amp[13] - amp[78] +
      Kokkos::complex<double> (0, 1) * amp[80] - Kokkos::complex<double> (0, 1)
      * amp[81] - Kokkos::complex<double> (0, 1) * amp[82] -
      Kokkos::complex<double> (0, 1) * amp[84] - amp[88] -
      Kokkos::complex<double> (0, 1) * amp[89] - Kokkos::complex<double> (0, 1)
      * amp[102] + Kokkos::complex<double> (0, 1) * amp[104] -
      Kokkos::complex<double> (0, 1) * amp[105] + amp[107] +
      Kokkos::complex<double> (0, 1) * amp[108] - amp[127] -
      Kokkos::complex<double> (0, 1) * amp[129] + Kokkos::complex<double> (0,
      1) * amp[131] + amp[132] - amp[134] + Kokkos::complex<double> (0, 1) *
      amp[139] + Kokkos::complex<double> (0, 1) * amp[140] + amp[145] +
      amp[146];
  jamp[19] = -amp[58] + Kokkos::complex<double> (0, 1) * amp[59] -
      Kokkos::complex<double> (0, 1) * amp[63] + Kokkos::complex<double> (0, 1)
      * amp[65] - Kokkos::complex<double> (0, 1) * amp[66] +
      Kokkos::complex<double> (0, 1) * amp[67] + amp[68] -
      Kokkos::complex<double> (0, 1) * amp[75] - amp[76] + amp[78] -
      Kokkos::complex<double> (0, 1) * amp[80] + Kokkos::complex<double> (0, 1)
      * amp[82] + Kokkos::complex<double> (0, 1) * amp[83] +
      Kokkos::complex<double> (0, 1) * amp[85] - amp[99] +
      Kokkos::complex<double> (0, 1) * amp[102] - Kokkos::complex<double> (0,
      1) * amp[104] + Kokkos::complex<double> (0, 1) * amp[105] - amp[107] -
      Kokkos::complex<double> (0, 1) * amp[108] + Kokkos::complex<double> (0,
      1) * amp[129] + Kokkos::complex<double> (0, 1) * amp[130] - amp[132] -
      amp[133] + amp[151] + amp[152] - Kokkos::complex<double> (0, 1) *
      amp[154] - Kokkos::complex<double> (0, 1) * amp[155];
  jamp[20] = +Kokkos::complex<double> (0, 1) * amp[1] + Kokkos::complex<double>
      (0, 1) * amp[2] + Kokkos::complex<double> (0, 1) * amp[3] +
      Kokkos::complex<double> (0, 1) * amp[4] - amp[10] -
      Kokkos::complex<double> (0, 1) * amp[11] + amp[13] - amp[61] -
      Kokkos::complex<double> (0, 1) * amp[62] + Kokkos::complex<double> (0, 1)
      * amp[63] + Kokkos::complex<double> (0, 1) * amp[64] +
      Kokkos::complex<double> (0, 1) * amp[66] - amp[68] +
      Kokkos::complex<double> (0, 1) * amp[69] + Kokkos::complex<double> (0, 1)
      * amp[111] - Kokkos::complex<double> (0, 1) * amp[113] +
      Kokkos::complex<double> (0, 1) * amp[114] + amp[118] -
      Kokkos::complex<double> (0, 1) * amp[119] - amp[128] -
      Kokkos::complex<double> (0, 1) * amp[130] - Kokkos::complex<double> (0,
      1) * amp[131] + amp[133] + amp[134] + Kokkos::complex<double> (0, 1) *
      amp[138] - Kokkos::complex<double> (0, 1) * amp[140] + amp[144] -
      amp[146];
  jamp[21] = -amp[39] + Kokkos::complex<double> (0, 1) * amp[41] - amp[42] -
      Kokkos::complex<double> (0, 1) * amp[43] - amp[48] + amp[52] + amp[53] +
      amp[61] + Kokkos::complex<double> (0, 1) * amp[62] -
      Kokkos::complex<double> (0, 1) * amp[63] - Kokkos::complex<double> (0, 1)
      * amp[64] - Kokkos::complex<double> (0, 1) * amp[66] + amp[68] -
      Kokkos::complex<double> (0, 1) * amp[69] + Kokkos::complex<double> (0, 1)
      * amp[102] + Kokkos::complex<double> (0, 1) * amp[103] +
      Kokkos::complex<double> (0, 1) * amp[105] + Kokkos::complex<double> (0,
      1) * amp[106] - amp[107] + Kokkos::complex<double> (0, 1) * amp[112] +
      Kokkos::complex<double> (0, 1) * amp[113] + Kokkos::complex<double> (0,
      1) * amp[115] + Kokkos::complex<double> (0, 1) * amp[129] +
      Kokkos::complex<double> (0, 1) * amp[130] - amp[132] - amp[133] +
      Kokkos::complex<double> (0, 1) * amp[157] + Kokkos::complex<double> (0,
      1) * amp[158];
  jamp[22] = +Kokkos::complex<double> (0, 1) * amp[0] + Kokkos::complex<double>
      (0, 1) * amp[1] + Kokkos::complex<double> (0, 1) * amp[3] +
      Kokkos::complex<double> (0, 1) * amp[5] + Kokkos::complex<double> (0, 1)
      * amp[9] - amp[10] - amp[17] + amp[58] - Kokkos::complex<double> (0, 1) *
      amp[59] + Kokkos::complex<double> (0, 1) * amp[63] -
      Kokkos::complex<double> (0, 1) * amp[65] + Kokkos::complex<double> (0, 1)
      * amp[66] - Kokkos::complex<double> (0, 1) * amp[67] - amp[68] - amp[98]
      - Kokkos::complex<double> (0, 1) * amp[100] + amp[101] +
      Kokkos::complex<double> (0, 1) * amp[120] - Kokkos::complex<double> (0,
      1) * amp[122] - Kokkos::complex<double> (0, 1) * amp[124] -
      Kokkos::complex<double> (0, 1) * amp[130] - Kokkos::complex<double> (0,
      1) * amp[131] + amp[133] + amp[134] + amp[150] - amp[152] -
      Kokkos::complex<double> (0, 1) * amp[153] + Kokkos::complex<double> (0,
      1) * amp[155];
  jamp[23] = -Kokkos::complex<double> (0, 1) * amp[0] - Kokkos::complex<double>
      (0, 1) * amp[1] - Kokkos::complex<double> (0, 1) * amp[3] -
      Kokkos::complex<double> (0, 1) * amp[5] - Kokkos::complex<double> (0, 1)
      * amp[9] + amp[10] + amp[17] + amp[39] - Kokkos::complex<double> (0, 1) *
      amp[41] - amp[47] - Kokkos::complex<double> (0, 1) * amp[49] + amp[50] +
      amp[51] - amp[53] - Kokkos::complex<double> (0, 1) * amp[102] -
      Kokkos::complex<double> (0, 1) * amp[103] - Kokkos::complex<double> (0,
      1) * amp[105] - Kokkos::complex<double> (0, 1) * amp[106] + amp[107] -
      Kokkos::complex<double> (0, 1) * amp[120] - Kokkos::complex<double> (0,
      1) * amp[121] - Kokkos::complex<double> (0, 1) * amp[123] -
      Kokkos::complex<double> (0, 1) * amp[129] + Kokkos::complex<double> (0,
      1) * amp[131] + amp[132] - amp[134] + Kokkos::complex<double> (0, 1) *
      amp[156] - Kokkos::complex<double> (0, 1) * amp[158];

  // Sum and square the color flows to get the matrix element
  // printf("1 ztemp = %7e,%7e\n",ztemp.real(),ztemp.imag());
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    // printf("2 ztemp = %7e,%7e\n",ztemp.real(),ztemp.imag());
    for(j = 0; j < ncolor; j++ ){
      // printf("3 ztemp = %7e,%7e cf[%d][%d] = %7e jamp[%d] = %7e,%7e matrix = %7e\n",
        // ztemp.real(),ztemp.imag(),i,j,cf[i][j],j,jamp[j].real(),jamp[j].imag(),matrix);
      ztemp = ztemp + cf[i][j] * jamp[j];
      // printf("4 ztemp = %7e,%7e cf[%d][%d] = %7e jamp[%d] = %7e,%7e matrix = %7e\n",
        // ztemp.real(),ztemp.imag(),i,j,cf[i][j],j,jamp[j].real(),jamp[j].imag(),matrix);
    }
    matrix = matrix + (ztemp * conj(jamp[i])).real()/denom[i]; 
  }

  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

} // end calculate_wavefunction()




//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
template <typename ExecSpace>
void sigmaKin(Kokkos::View<double***,ExecSpace> momenta, 
    Kokkos::View<double*,ExecSpace> output,
    Kokkos::View<int**,ExecSpace> cHel,
    Kokkos::View<double*,ExecSpace> cIPD,
    Kokkos::View<Kokkos::complex<double>*,ExecSpace> cIPC,
    const int range) 
{
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();

  // Reset color flows

  Kokkos::View<Kokkos::complex<double>*,ExecSpace> amp("amp",159);

  // for (int xx = 0; xx < 384; ++xx) {

  Kokkos::parallel_for("sigmaKin",Kokkos::RangePolicy<ExecSpace>(0,range),
    KOKKOS_LAMBDA(const int& tid){

      const int nprocesses = 1;

      // char *devPtr = (char *)tp.ptr;
      // size_t dpt = tp.pitch;
      // size_t slicePitch = dpt * 6;

      // char *dps = devPtr + dim * slicePitch;
      double matrix_element[nprocesses]; 

      // double local_m[6][3];

      // for (int i=0; i<20;i++){
      // printf(" %f ", allmomenta[i]);
      // }
      // printf("\n");
      // printf("DIM is %i/%i\n", tid, DIM);
      // for (int i = 0; i < 6; i++ )
      // {
      //   for (int j = 0; j < 3; j++ )
      //   {
      //     local_m[i][j] = momenta(tid,i,j);
      //     // printf(" %f ", local_m[i][j]);
      //   }
      //   // printf("\n");
      // }


      // Local variables and constants
      const int ncomb = 64;
      // static bool goodhel[ncomb] = {ncomb * false};
      // static int ntry = 0, sum_hel = 0, ngood = 0;
      // static int igood[ncomb];
      // static int jhel;
      // std::complex<double> **wfs;
      // double t[1];
      // Helicities for the process
      // static const int helicities[ncomb][nexternal] =
      // {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,1},{-1,-1,-1,-1,1,-1},{-1,-1,-1,-1,1,1
      // },{-1,-1,-1,1,-1,-1},{-1,-1,-1,1,-1,1},{-1,-1,-1,1,1,-1},{-1,-1,-1,1,1,1},{
      // -1,-1,1,-1,-1,-1},{-1,-1,1,-1,-1,1},{-1,-1,1,-1,1,-1},{-1,-1,1,-1,1,1},{-1,
      // -1,1,1,-1,-1},{-1,-1,1,1,-1,1},{-1,-1,1,1,1,-1},{-1,-1,1,1,1,1},{-1,1,-1,-1
      // ,-1,-1},{-1,1,-1,-1,-1,1},{-1,1,-1,-1,1,-1},{-1,1,-1,-1,1,1},{-1,1,-1,1,-1,
      // -1},{-1,1,-1,1,-1,1},{-1,1,-1,1,1,-1},{-1,1,-1,1,1,1},{-1,1,1,-1,-1,-1},{-1
      // ,1,1,-1,-1,1},{-1,1,1,-1,1,-1},{-1,1,1,-1,1,1},{-1,1,1,1,-1,-1},{-1,1,1,1,-
      // 1,1},{-1,1,1,1,1,-1},{-1,1,1,1,1,1},{1,-1,-1,-1,-1,-1},{1,-1,-1,-1,-1,1},{1
      // ,-1,-1,-1,1,-1},{1,-1,-1,-1,1,1},{1,-1,-1,1,-1,-1},{1,-1,-1,1,-1,1},{1,-1,-
      // 1,1,1,-1},{1,-1,-1,1,1,1},{1,-1,1,-1,-1,-1},{1,-1,1,-1,-1,1},{1,-1,1,-1,1,-
      // 1},{1,-1,1,-1,1,1},{1,-1,1,1,-1,-1},{1,-1,1,1,-1,1},{1,-1,1,1,1,-1},{1,-1,1
      // ,1,1,1},{1,1,-1,-1,-1,-1},{1,1,-1,-1,-1,1},{1,1,-1,-1,1,-1},{1,1,-1,-1,1,1}
      // ,{1,1,-1,1,-1,-1},{1,1,-1,1,-1,1},{1,1,-1,1,1,-1},{1,1,-1,1,1,1},{1,1,1,-1,
      // -1,-1},{1,1,1,-1,-1,1},{1,1,1,-1,1,-1},{1,1,1,-1,1,1},{1,1,1,1,-1,-1},{1,1,
      // 1,1,-1,1},{1,1,1,1,1,-1},{1,1,1,1,1,1}};
      // Denominators: spins, colors and identical particles
      const int denominators[1] = {512};


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

      auto local_mom = Kokkos::subview(momenta,tid,Kokkos::ALL,Kokkos::ALL);
      for (int ihel = 0; ihel < ncomb; ++ihel)
      {
        auto local_cHel = Kokkos::subview(cHel,ihel,Kokkos::ALL);
        // printf("tid = %d ihel = %d cHel[%d][0] = %d matrix = %7e\n",tid,ihel,ihel,local_cHel(0),matrix_element[0]);
        calculate_wavefunctions(local_cHel, local_mom, cIPD, cIPC, matrix_element[0]);
        // printf("tid = %d ihel = %d matrix = %7e\n",tid,ihel,matrix_element[0]);
        
      }

      // printf("tid = %d  me = %f\n",tid,matrix_element[0]);

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
    });

}

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ g g WEIGHTED<=4 @1
//--------------------------------------------------------------------------

template <class ExecSpace>
class CPPProcess
{
  public:

    CPPProcess(int numiterations, int gpublocks, int gputhreads,
        bool verbose = false, bool debug = false): m_numiterations(numiterations),
        gpu_nblocks(gpublocks), gpu_nthreads(gputhreads), dim(gpu_nblocks * gpu_nthreads),
        cHel("cHel",ncomb,nexternal), hHel("hHel",ncomb,nexternal), 
        cmME("cmME",nexternal), hmME("hmME",nexternal),
        cIPC("cIPC",nexternal), hIPC("hIPC",nexternal), 
        cIPD("cIPD",ninitial), hIPD("hIPD",ninitial) {

      const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1}, {-1, -1,
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
        {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1},
      };
      for(int i=0;i<ncomb;++i)
        for(int j=0;j<nexternal;++j){
            hHel(i,j) = tHel[i][j];
        }
      Kokkos::deep_copy(cHel,hHel);
    }

    ~CPPProcess() = default;

    //--------------------------------------------------------------------------
    // Initialize process.
    void initProc(const std::string param_card_name)
    {
        // Helicities for the process - nodim
        // Kokkos::parallel_for("CPPProcessInit",Kokkos::RangePolicy<ExecSpace>(0,ncomb),
        // KOKKOS_CLASS_LAMBDA(int i){
        //   for(auto j = 0; j < nexternal; ++j)
        //     cHel(i,j) = (i >> j) & 1;
        // });

        // Instantiate the model class and set parameters that stay fixed during run
        pars = Parameters_sm::getInstance(); 
        SLHAReader slha(param_card_name); 
        pars->setIndependentParameters(slha); 
        pars->setIndependentCouplings(); 
        pars->printIndependentParameters(); 
        pars->printIndependentCouplings(); 
        pars->setDependentParameters(); 
        pars->setDependentCouplings(); 
        // Set external particle masses for this matrix element
        hmME(0) = pars->ZERO;
        hmME(1) = pars->ZERO;
        hmME(2) = pars->mdl_MT;
        hmME(3) = pars->mdl_MT;
        hmME(4) = pars->ZERO;
        hmME(5) = pars->ZERO;
        Kokkos::deep_copy(cmME,hmME);

        hIPC(0) = pars->GC_10;
        hIPC(1) = pars->GC_11;
        hIPC(2) = pars->GC_12;
        Kokkos::deep_copy(cIPC,hIPC);

        hIPD(0) = pars->mdl_MT;
        hIPD(1) = pars->mdl_WT; 
        Kokkos::deep_copy(cIPD,hIPD);
    }

    // virtual int code() const {return 1;}

    // void setInitial(int inid1, int inid2)
    // {
    //   id1 = inid1;
    //   id2 = inid2;
    // }

    // int getDim() const {return dim;}

    // int getNIOParticles() const {return nexternal;}


    // Constants for array limits
    static const int ninitial = 2; 
    static const int nprocesses = 1; 
    static const int ncomb = 64; 
    static const int nexternal = 6; 

    Kokkos::View<int**,ExecSpace> cHel;
    typename Kokkos::View<int**,ExecSpace>::HostMirror hHel;

    Kokkos::View<double*,ExecSpace> cmME;
    typename Kokkos::View<double*,ExecSpace>::HostMirror hmME;

    Kokkos::View<Kokkos::complex<double>*,ExecSpace> cIPC;
    typename Kokkos::View<Kokkos::complex<double>*,ExecSpace>::HostMirror hIPC;

    Kokkos::View<double*,ExecSpace> cIPD;
    typename Kokkos::View<double*,ExecSpace>::HostMirror hIPD;

  private:
    int m_numiterations; 
    // gpu variables
    int gpu_nblocks; 
    int gpu_nthreads; 
    int dim;  // gpu_nblocks * gpu_nthreads;
    // vector with external particle masses
    // std::vector<double> mME;
    // Pointer to the model parameters
    Parameters_sm * pars; 

    // // print verbose info
    // bool m_verbose; 

    // // print debug info
    // bool m_debug; 

    // static const int nwavefuncs = 6; 
    // static const int namplitudes = 159; 
    // static const int wrows = 63; 
    // // static const int nioparticles = 6;

    // Kokkos::complex<double> * * amp; 




    // // Initial particle ids
    // int id1, id2; 

};


#endif  // MG5_Sigma_sm_gg_ttxgg_H
