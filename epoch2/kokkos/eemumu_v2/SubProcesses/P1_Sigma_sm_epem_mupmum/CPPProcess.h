//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#include "Parameters_sm.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_Complex.hpp"
#include <ctime>


template <typename Layout, typename Device, typename MemSpace, typename ExecSpace>
void ixxxxx(Kokkos::View<double**,Layout,Device,MemSpace> momenta,
            double fmass,
            const int nhel, 
            const int nsf, 
            Kokkos::View<Kokkos::complex<double>**,Layout,Device,MemSpace> w, // fi[6]
            const int num_events)
{
  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){

    Kokkos::complex<double> chi[2]; 
    double sf[2], sfomega[2], omega[2], pp = 0., pp3 = 0., sqp0p3 = 0., sqm[2]; 
    int ip = 0, im = 0, nh = 0;

    auto p = Kokkos::subview(momenta,tid,Kokkos::ALL);
    auto fi = Kokkos::subview(w,tid,Kokkos::ALL);

    double p2 = p(1) * p(1) + p(2) * p(2) + p(3) * p(3);
    double p_sqrt = sqrt(p2);

    fi(0) = Kokkos::complex<double> (-p(0) * nsf, -p(3) * nsf); 
    fi(1) = Kokkos::complex<double> (-p(1) * nsf, -p(2) * nsf); 
    nh = nhel * nsf; 
    if (fmass != 0.0)
    {
      pp = min(p[0], p_sqrt);
      if (pp == 0.0)
      {
        sqm[0] = sqrt(std::abs(fmass)); 
        sqm[1] = (fmass < 0) ? - abs(sqm[0]) : abs(sqm[0]); 
        ip = (1 + nh)/2; 
        im = (1 - nh)/2; 
        fi(2) = ip * sqm[ip]; 
        fi(3) = im * nsf * sqm[ip]; 
        fi(4) = ip * nsf * sqm[im]; 
        fi(5) = im * sqm[im]; 
      }
      else
      {
        sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
        sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
        omega[0] = sqrt(p(0) + pp); 
        omega[1] = fmass/omega[0]; 
        ip = (1 + nh)/2; 
        im = (1 - nh)/2; 
        sfomega[0] = sf[0] * omega[ip]; 
        sfomega[1] = sf[1] * omega[im]; 
        pp3 = max(pp + p(3), 0.0); 
        chi[0] = Kokkos::complex<double> (sqrt(pp3 * 0.5/pp), 0); 
        if (pp3 == 0.0)
        {
          chi[1] = Kokkos::complex<double> (-nh, 0); 
        }
        else
        {
          chi[1] = 
          Kokkos::complex<double> (nh * p(1), p(2))/sqrt(2.0 * pp * pp3); 
        }
        fi(2) = sfomega[0] * chi[im]; 
        fi(3) = sfomega[0] * chi[ip]; 
        fi(4) = sfomega[1] * chi[im]; 
        fi(5) = sfomega[1] * chi[ip]; 
      }
    }
    else
    {
      if (p(1) == 0.0 and p(2) == 0.0 and p(3) < 0.0)
      {
        sqp0p3 = 0.0; 
      }
      else
      {
        sqp0p3 = sqrt(max(p(0) + p(3), 0.0)) * nsf; 
      }
      chi[0] = Kokkos::complex<double> (sqp0p3, 0.0); 
      if (sqp0p3 == 0.0)
      {
        chi[1] = Kokkos::complex<double> (-nhel * sqrt(2.0 * p(0)), 0.0); 
      }
      else
      {
        chi[1] = Kokkos::complex<double> (nh * p(1), p(2))/sqp0p3; 
      }
      if (nh == 1)
      {
        fi(2) = Kokkos::complex<double> (0.0, 0.0); 
        fi(3) = Kokkos::complex<double> (0.0, 0.0); 
        fi(4) = chi[0]; 
        fi(5) = chi[1]; 
      }
      else
      {
        fi(2) = chi[1]; 
        fi(3) = chi[0]; 
        fi(4) = Kokkos::complex<double> (0.0, 0.0); 
        fi(5) = Kokkos::complex<double> (0.0, 0.0); 
      }
    }
  }); // end parallel_for
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
  // double p[4] = {0, pvec(1), pvec(2), pvec(3)};
  // p[0] = sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]+fmass*fmass);
  double p[4] = {0, 0, 0, 0}; 
  printf("scalar not supported so far. to do: fix mass issue"); 
  sc[2] = Kokkos::complex<double> (1.00, 0.00); 
  sc[0] = Kokkos::complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = Kokkos::complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}


template <typename Layout, typename Device, typename MemSpace, typename ExecSpace>
void oxxxxx(Kokkos::View<double**,Layout,Device,MemSpace> momenta,
            double fmass,
            const int nhel, 
            const int nsf, 
            Kokkos::View<Kokkos::complex<double>**,Layout,Device,MemSpace> w, // fo[6]
            const int num_events) 
{
  // Kokkos::Timer total_timer;
  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){
    Kokkos::complex<double> chi[2];
    double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
    int nh, ip, im;

    auto p = Kokkos::subview(momenta,tid,Kokkos::ALL);
    auto fo = Kokkos::subview(w,tid,Kokkos::ALL);

    double p2 = p(1) * p(1) + p(2) * p(2) + p(3) * p(3);
    double p_sqrt = sqrt(p2);
    double pf = p2 + fmass * fmass;

    fo(0) = Kokkos::complex<double> (pf * nsf, p(3) * nsf); 
    fo(1) = Kokkos::complex<double> (p(1) * nsf, p(2) * nsf); 
    nh = nhel * nsf;
    if (fmass != 0.000)
    {
      pp = min(pf, p_sqrt); 
      if (pp == 0.000)
      {
        sqm[0] = sqrt(std::abs(fmass)); 
        sqm[1] = (fmass < 0) ? - abs(sqm[0]) : abs(sqm[0]); 
        ip = -((1 - nh)/2) * nhel; 
        im = (1 + nh)/2 * nhel; 
        fo(2) = im * sqm[std::abs(ip)]; 
        fo(3) = ip * nsf * sqm[std::abs(ip)]; 
        fo(4) = im * nsf * sqm[std::abs(im)]; 
        fo(5) = ip * sqm[std::abs(im)];
      }
      else
      {
        pp = min(pf, p_sqrt); 
        sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
        sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
        omega[0] = sqrt(pf + pp); 
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
          Kokkos::complex<double> (nh * p(1), -p(2))/sqrt(2.0 * pp * pp3); 
        }
        fo(2) = sfomeg[1] * chi[im]; 
        fo(3) = sfomeg[1] * chi[ip]; 
        fo(4) = sfomeg[0] * chi[im]; 
        fo(5) = sfomeg[0] * chi[ip]; 
      }
    }
    else
    {
      if ((p(1) == 0.00) and (p(2) == 0.00) and (p(3) < 0.00))
      {
        sqp0p3 = 0.00; 
      }
      else
      {
        sqp0p3 = sqrt(max(pf + p(3), 0.00)) * nsf; 
      }
      chi[0] = Kokkos::complex<double> (sqp0p3, 0.00); 
      if (sqp0p3 == 0.000)
      {
        chi[1] = Kokkos::complex<double> (-nhel, 0.00) * sqrt(2.0 * pf); 
      }
      else
      {
        chi[1] = Kokkos::complex<double> (nh * p(1), -p(2))/sqp0p3; 
      }
      if (nh == 1)
      {
        fo(2) = chi[0]; 
        fo(3) = chi[1]; 
        fo(4) = Kokkos::complex<double> (0.00, 0.00); 
        fo(5) = Kokkos::complex<double> (0.00, 0.00); 
      }
      else
      {
        fo(2) = Kokkos::complex<double> (0.00, 0.00); 
        fo(3) = Kokkos::complex<double> (0.00, 0.00); 
        fo(4) = chi[1]; 
        fo(5) = chi[0]; 
      }
    }
  }); // end parallel_for
  // ExecSpace().fence();
  // printf("oxxxx total: %e\n",total_timer.seconds());
  return;
}


KOKKOS_FUNCTION void FFV2_0(const Kokkos::complex<double> F1[], const
Kokkos::complex<double> F2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> COUP, Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  Kokkos::complex<double> TMP0; 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))); 
  (*vertex) = COUP * - cI * TMP0; 
}


KOKKOS_FUNCTION void FFV2_3(const Kokkos::complex<double> F1[], const
Kokkos::complex<double> F2[], const Kokkos::complex<double> COUP, const
double M3, const double W3, Kokkos::complex<double> V3[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> denom; 
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


KOKKOS_FUNCTION void FFV4_0(const Kokkos::complex<double> F1[], const
Kokkos::complex<double> F2[], const Kokkos::complex<double> V3[], const
Kokkos::complex<double> COUP, Kokkos::complex<double> * vertex)
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  Kokkos::complex<double> TMP0; 
  Kokkos::complex<double> TMP2; 
  TMP2 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) + 
  F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5]))); 
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) + 
  F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5]))); 
  (*vertex) = COUP * (-1.) * (+cI * (TMP0) + 2. * cI * (TMP2)); 
}


KOKKOS_FUNCTION void FFV4_3(const Kokkos::complex<double> F1[], const
Kokkos::complex<double> F2[], const Kokkos::complex<double> COUP, const
double M3, const double W3, Kokkos::complex<double> V3[])
{
  Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
  double OM3; 
  double P3[4]; 
  Kokkos::complex<double> TMP1; 
  Kokkos::complex<double> TMP3; 
  Kokkos::complex<double> denom; 
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


template <typename Stride, typename Device, typename MemSpace, typename ExecSpace>
void FFV1_0(Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_F1,
            Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_F2, 
            Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_V3, 
            const Kokkos::complex<double> COUP,
            Kokkos::View<Kokkos::complex<double>*,Stride,Device,MemSpace> amp_vertex,
            const int num_events)
{
  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){
    auto F1 = Kokkos::subview(w_F1,tid,Kokkos::ALL);
    auto F2 = Kokkos::subview(w_F2,tid,Kokkos::ALL);
    auto V3 = Kokkos::subview(w_V3,tid,Kokkos::ALL);
    auto vertex = Kokkos::subview(amp_vertex,tid);
    Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
    Kokkos::complex<double> TMP4;
    TMP4 = (F1(2) * (F2(4) * (V3(2) + V3(5)) + F2(5) * (V3(3) + cI * (V3(4)))) + 
    (F1(3) * (F2(4) * (V3(3) - cI * (V3(4))) + F2(5) * (V3(2) - V3(5))) + 
    (F1(4) * (F2(2) * (V3(2) - V3(5)) - F2(3) * (V3(3) + cI * (V3(4)))) + 
    F1(5) * (F2(2) * (-V3(3) + cI * (V3(4))) + F2(3) * (V3(2) + V3(5)))))); 
    vertex() = COUP * - cI * TMP4; 
  });
}


template <typename Layout, typename Device, typename MemSpace, typename ExecSpace>
void FFV1P0_3(Kokkos::View<Kokkos::complex<double>**,Layout,Device,MemSpace> w_F1,
              Kokkos::View<Kokkos::complex<double>**,Layout,Device,MemSpace> w_F2,
              const Kokkos::complex<double> COUP, 
              const double M3,
              const double W3, 
              Kokkos::View<Kokkos::complex<double>**,Layout,Device,MemSpace> w_V3,
              const int num_events)
{

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){

    auto F1 = Kokkos::subview(w_F1,tid,Kokkos::ALL);
    auto F2 = Kokkos::subview(w_F2,tid,Kokkos::ALL);
    auto V3 = Kokkos::subview(w_V3,tid,Kokkos::ALL);

    Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.);
    double P3[4];
    Kokkos::complex<double> denom;
    V3(0) = +F1(0) + F2(0); 
    V3(1) = +F1(1) + F2(1); 
    P3[0] = -V3(0).real(); 
    P3[1] = -V3(1).real(); 
    P3[2] = -V3(1).imag(); 
    P3[3] = -V3(0).imag(); 
    denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * 
    P3[3]) - M3 * (M3 - cI * W3)); 
    V3(2) = denom * (-cI) * (F1(2) * F2(4) + F1(3) * F2(5) + F1(4) * F2(2) + 
    F1(5) * F2(3)); 
    V3(3) = denom * (-cI) * (-F1(2) * F2(5) - F1(3) * F2(4) + F1(4) * F2(3) + 
    F1(5) * F2(2)); 
    V3(4) = denom * (-cI) * (-cI * (F1(2) * F2(5) + F1(5) * F2(2)) + cI * (F1(3)
     * F2(4) + F1(4) * F2(3))); 
    V3(5) = denom * (-cI) * (-F1(2) * F2(4) - F1(5) * F2(3) + F1(3) * F2(5) + 
    F1(4) * F2(2)); 
  });
}


template <typename Stride, typename Device, typename MemSpace, typename ExecSpace>
void FFV2_4_0(Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_F1,
              Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_F2, 
              Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_V3, 
              const Kokkos::complex<double> COUP1, 
              const Kokkos::complex<double> COUP2, 
              Kokkos::View<Kokkos::complex<double>*,Stride,Device,MemSpace> amp_vertex,
              const int num_events)
{
  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){
    auto F1 = Kokkos::subview(w_F1,tid,Kokkos::ALL);
    auto F2 = Kokkos::subview(w_F2,tid,Kokkos::ALL);
    auto V3 = Kokkos::subview(w_V3,tid,Kokkos::ALL);
    auto vertex = Kokkos::subview(amp_vertex,tid);

    Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
    Kokkos::complex<double> TMP0; 
    Kokkos::complex<double> TMP2; 
    TMP2 = (F1(4) * (F2(2) * (V3(2) - V3(5)) - F2(3) * (V3(3) + cI * (V3(4)))) + 
    F1(5) * (F2(2) * (-V3(3) + cI * (V3(4))) + F2(3) * (V3(2) + V3(5)))); 
    TMP0 = (F1(2) * (F2(4) * (V3(2) + V3(5)) + F2(5) * (V3(3) + cI * (V3(4)))) + 
    F1(3) * (F2(4) * (V3(3) - cI * (V3(4))) + F2(5) * (V3(2) - V3(5)))); 
    vertex() = (-1.) * (COUP2 * (+cI * (TMP0) + 2. * cI * (TMP2)) + cI * (TMP0 * 
    COUP1)); 
  });
}


template <typename Stride, typename Device, typename MemSpace, typename ExecSpace>
void FFV2_4_3(Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_F1,
              Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_F2,
              const Kokkos::complex<double> COUP1, 
              const Kokkos::complex<double> COUP2, 
              const double M3, 
              const double W3, 
              Kokkos::View<Kokkos::complex<double>**,Stride,Device,MemSpace> w_V3,
              const int num_events)
{

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){

    auto F1 = Kokkos::subview(w_F1,tid,Kokkos::ALL);
    auto F2 = Kokkos::subview(w_F2,tid,Kokkos::ALL);
    auto V3 = Kokkos::subview(w_V3,tid,Kokkos::ALL);

    Kokkos::complex<double> cI = Kokkos::complex<double> (0., 1.); 
    double OM3; 
    double P3[4]; 
    Kokkos::complex<double> TMP1; 
    Kokkos::complex<double> TMP3; 
    Kokkos::complex<double> denom; 
    OM3 = 0.; 
    if (M3 != 0.)
      OM3 = 1./(M3 * M3); 
    V3(0) = +F1(0) + F2(0); 
    V3(1) = +F1(1) + F2(1); 
    P3[0] = -V3(0).real(); 
    P3[1] = -V3(1).real(); 
    P3[2] = -V3(1).imag(); 
    P3[3] = -V3(0).imag(); 
    TMP1 = (F1(2) * (F2(4) * (P3[0] + P3[3]) + F2(5) * (P3[1] + cI * (P3[2]))) + 
    F1(3) * (F2(4) * (P3[1] - cI * (P3[2])) + F2(5) * (P3[0] - P3[3]))); 
    TMP3 = (F1(4) * (F2(2) * (P3[0] - P3[3]) - F2(3) * (P3[1] + cI * (P3[2]))) + 
    F1(5) * (F2(2) * (-P3[1] + cI * (P3[2])) + F2(3) * (P3[0] + P3[3]))); 
    denom = 1./((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] * 
    P3[3]) - M3 * (M3 - cI * W3)); 
    V3(2) = denom * (-2. * cI) * (COUP2 * (OM3 * - 1./2. * P3[0] * (TMP1 + 2. * 
    (TMP3)) + (+1./2. * (F1(2) * F2(4) + F1(3) * F2(5)) + F1(4) * F2(2) + 
    F1(5) * F2(3))) + 1./2. * (COUP1 * (F1(2) * F2(4) + F1(3) * F2(5) - P3[0]
     * OM3 * TMP1))); 
    V3(3) = denom * (-2. * cI) * (COUP2 * (OM3 * - 1./2. * P3[1] * (TMP1 + 2. * 
    (TMP3)) + (-1./2. * (F1(2) * F2(5) + F1(3) * F2(4)) + F1(4) * F2(3) + 
    F1(5) * F2(2))) - 1./2. * (COUP1 * (F1(2) * F2(5) + F1(3) * F2(4) + P3[1]
     * OM3 * TMP1))); 
    V3(4) = denom * cI * (COUP2 * (OM3 * P3[2] * (TMP1 + 2. * (TMP3)) + (+cI * 
    (F1(2) * F2(5)) - cI * (F1(3) * F2(4)) - 2. * cI * (F1(4) * F2(3)) + 2. * 
    cI * (F1(5) * F2(2)))) + COUP1 * (+cI * (F1(2) * F2(5)) - cI * (F1(3) * 
    F2(4)) + P3[2] * OM3 * TMP1)); 
    V3(5) = denom * 2. * cI * (COUP2 * (OM3 * 1./2. * P3[3] * (TMP1 + 2. * 
    (TMP3)) + (+1./2. * (F1(2) * F2(4)) - 1./2. * (F1(3) * F2(5)) - F1(4) * 
    F2(2) + F1(5) * F2(3))) + 1./2. * (COUP1 * (F1(2) * F2(4) + P3[3] * OM3 * 
    TMP1 - F1(3) * F2(5)))); 
  });
}

//       Kokkos::View<int *, Kokkos::LayoutStride, Kokkos::Device<Kokkos::CudaUVMSpace::execution_space, Kokkos::CudaUVMSpace::memory_space>, Kokkos::MemoryManaged>,
// const Kokkos::View<double ***, Kokkos::DefaultExecutionSpace>, 
// const Kokkos::View<double *, Kokkos::DefaultExecutionSpace>, 
// const Kokkos::View<Kokkos::complex<double> *, Kokkos::DefaultExecutionSpace>, 
// Kokkos::View<double *, Kokkos::Cuda::array_layout, Kokkos::Device<Kokkos::CudaUVMSpace::execution_space, Kokkos::CudaUVMSpace::memory_space>, Kokkos::MemoryManaged>, 
// Kokkos::View<Kokkos::complex<double> **, Kokkos::DefaultExecutionSpace>, 
// Kokkos::View<Kokkos::complex<double> ***, Kokkos::DefaultExecutionSpace>, 
// const int)
// Evaluate |M|^2 for each subprocess
template <typename Stride, typename Layout, typename Device, typename MemSpace, typename ExecSpace>
void calculate_wavefunctions(
    const Kokkos::View<int*,Stride,Device,MemSpace>& cHel,
    const Kokkos::View<double***,ExecSpace>& momenta,
    const Kokkos::View<double*,ExecSpace>& cIPD,
    const Kokkos::View<Kokkos::complex<double>*,ExecSpace>& cIPC,
    Kokkos::View<double*,Layout,Device,MemSpace>& matrix,
    Kokkos::View<Kokkos::complex<double>**,ExecSpace>& amp,
    Kokkos::View<Kokkos::complex<double>***,ExecSpace>& w,
    const int& num_events)
{
  // Kokkos::Timer total_timer;
  // Kokkos::Timer timer;

  oxxxxx<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(momenta,Kokkos::ALL,0,Kokkos::ALL), 0., cHel(0), -1, 
         Kokkos::subview(w,Kokkos::ALL,0,Kokkos::ALL),num_events); 
  // ExecSpace().fence();
  // printf("oxA: %e ",timer.seconds());
  // timer.reset();
  ixxxxx<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(momenta,Kokkos::ALL,1,Kokkos::ALL), 0., cHel(1), -1,
         Kokkos::subview(w,Kokkos::ALL,1,Kokkos::ALL),num_events);  
  // ExecSpace().fence();
  // printf("ixA: %e ",timer.seconds());
  // timer.reset();
  ixxxxx<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(momenta,Kokkos::ALL,2,Kokkos::ALL), 0., cHel(2), -1,
         Kokkos::subview(w,Kokkos::ALL,2,Kokkos::ALL),num_events);  
  // ExecSpace().fence();
  // printf("ixB: %e ",timer.seconds());
  // timer.reset();
  oxxxxx<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(momenta,Kokkos::ALL,3,Kokkos::ALL), 0., cHel(3), -1,
         Kokkos::subview(w,Kokkos::ALL,3,Kokkos::ALL),num_events); 
  ExecSpace().fence();
  // printf("oxB: %e ",timer.seconds());
  // timer.reset();
  

  // Amplitude(s) for diagram number 1
  FFV1P0_3<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(w,Kokkos::ALL,1,Kokkos::ALL), 
           Kokkos::subview(w,Kokkos::ALL,0,Kokkos::ALL), cIPC(0), 0., 0., 
           Kokkos::subview(w,Kokkos::ALL,4,Kokkos::ALL),num_events);
  ExecSpace().fence();
  // printf("F1_A: %e ",timer.seconds());
  // timer.reset();
  FFV1_0<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(w,Kokkos::ALL,2,Kokkos::ALL), 
         Kokkos::subview(w,Kokkos::ALL,3,Kokkos::ALL), 
         Kokkos::subview(w,Kokkos::ALL,4,Kokkos::ALL), cIPC(0), 
         Kokkos::subview(amp,Kokkos::ALL,0),num_events);
  ExecSpace().fence();
  // Amplitude(s) for diagram number 2
  // printf("F1_B: %e ",timer.seconds());
  // timer.reset();
  FFV2_4_3<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(w,Kokkos::ALL,1,Kokkos::ALL), 
           Kokkos::subview(w,Kokkos::ALL,0,Kokkos::ALL), cIPC(1), cIPC(2), cIPD(0), cIPD(1), 
           Kokkos::subview(w,Kokkos::ALL,4,Kokkos::ALL),num_events);
  ExecSpace().fence();
  // printf("F2_A: %e ",timer.seconds());
  // timer.reset();
  FFV2_4_0<Layout,Device,MemSpace,ExecSpace>(Kokkos::subview(w,Kokkos::ALL,2,Kokkos::ALL), 
           Kokkos::subview(w,Kokkos::ALL,3,Kokkos::ALL), 
           Kokkos::subview(w,Kokkos::ALL,4,Kokkos::ALL), cIPC(1), cIPC(2),
           Kokkos::subview(amp,Kokkos::ALL,1),num_events);
  ExecSpace().fence();
  // printf("F2_B: %e \n",timer.seconds());
  // timer.reset();
  // ExecSpace().fence();

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){
    // const int ngraphs = 2;
    const int ncolor = 1; 
    Kokkos::complex<double> ztemp; 
    Kokkos::complex<double> jamp[ncolor]; 
    // The color matrix;
    static const double denom[ncolor] = {1}; 
    static const double cf[ncolor][ncolor] = {{1}}; 

    // Calculate color flows
    jamp[0] = -amp(tid,0) - amp(tid,1);

    // Sum and square the color flows to get the matrix element
    for(int i = 0; i < ncolor; i++ )
    {
      ztemp = 0.; 
      for(int j = 0; j < ncolor; j++ )
        ztemp = ztemp + cf[i][j] * jamp[j]; 
      matrix[tid] = matrix[tid] + (ztemp * conj(jamp[i])).real()/denom[i]; 
    }
  });
  ExecSpace().fence();
  // printf("sigmaKin loop: %e \n",timer.seconds());
  // printf("total calcw: %e \n",total_timer.seconds());
  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

}


//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
template <typename ExecSpace>
void sigmaKin(const Kokkos::View<double***,ExecSpace>& momenta, 
    Kokkos::View<double*,ExecSpace>& output,
    const Kokkos::View<int**,ExecSpace>& cHel,
    const Kokkos::View<double*,ExecSpace>& cIPD,
    const Kokkos::View<Kokkos::complex<double>*,ExecSpace>& cIPC,
    Kokkos::View<Kokkos::complex<double>**,ExecSpace>& amp,
    Kokkos::View<Kokkos::complex<double>***,ExecSpace>& w,
    Kokkos::View<double**,ExecSpace>& matrix_element,
    const int& num_events)
{
  // Kokkos::Timer timer,total_timer;
  
  const int ncomb = 16;

  // float subview = 0.;
  // float calcw = 0;

  for (int ihel = 0; ihel < ncomb; ++ihel)
  {
    // timer.reset();
    const auto local_cHel = Kokkos::subview(cHel,ihel,Kokkos::ALL);
    auto local_me = Kokkos::subview(matrix_element,Kokkos::ALL,0);
    // subview += timer.seconds();
    // printf("tid = %d ihel = %d cHel[%d][0] = %d matrix = %7e\n",tid,ihel,ihel,local_cHel(0),matrix_element[0]);
    // timer.reset();
    calculate_wavefunctions(local_cHel, momenta, cIPD, cIPC, 
                            local_me,
                            amp, w, num_events);
    // calcw += timer.seconds();
    // printf("tid = %d ihel = %d matrix = %7e\n",tid,ihel,matrix_element[0]);
    
  }

  // printf("subview cHel takes total %f per ihel %f \n",subview,subview/ncomb);
  // printf("calcw takes total %f per ihel %f \n",calcw,calcw/ncomb);
  // timer.reset();
  Kokkos::parallel_for("sigmaKin",Kokkos::RangePolicy<ExecSpace>(0,num_events),
  KOKKOS_LAMBDA(const int& tid){
    output[tid] = matrix_element(tid,0) / 4;
  }); // end parallel for

  ExecSpace().fence();

  // printf("sigmaKin loop: %f\n",timer.seconds());
  // printf("sigmaKin function time: %f\n",total_timer.seconds());
}

//==========================================================================
// A class for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1
//--------------------------------------------------------------------------
template <class ExecSpace>
class CPPProcess
{
  public:

    CPPProcess(int numiterations, int events_per_iter): 
        m_numiterations(numiterations), events_per_iter(events_per_iter),
        cHel("cHel",ncomb,nexternal), hHel("hHel",ncomb,nexternal), 
        cmME("cmME",nexternal), hmME("hmME",nexternal),
        cIPC("cIPC",nexternal), hIPC("hIPC",nexternal), 
        cIPD("cIPD",ninitial), hIPD("hIPD",ninitial) {


      // Helicities for the process - nodim
      static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1, -1, 1},
          {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1}, {-1, 1,
          1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1, 1, -1},
          {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1, 1, 1,
          1}};

      for(int i=0;i<ncomb;++i)
        for(int j=0;j<nexternal;++j){
            hHel(i,j) = tHel[i][j];
        }
      Kokkos::deep_copy(cHel,hHel);
    }

    ~CPPProcess() = default; 

    // Initialize process.
    void initProc(std::string param_card_name) 
    {
      // Instantiate the model class and set parameters that stay fixed during run
      pars = Parameters_sm::getInstance(); 
      SLHAReader slha(param_card_name); 
      pars->setIndependentParameters(slha); 
      pars->setIndependentCouplings(); 
      pars->setDependentParameters(); 
      pars->setDependentCouplings(); 
      // Set external particle masses for this matrix element
      hmME(0) = pars->ZERO; 
      hmME(1) = pars->ZERO; 
      hmME(2) = pars->ZERO; 
      hmME(3) = pars->ZERO;
      Kokkos::deep_copy(cmME,hmME);

      hIPC(0) = pars->GC_3;
      hIPC(1) = pars->GC_50;
      hIPC(2) = pars->GC_59;
      Kokkos::deep_copy(cIPC,hIPC);

      hIPD(0) = pars->mdl_MZ;
      hIPD(1) = pars->mdl_WZ;
      Kokkos::deep_copy(cIPD,hIPD);
    } 


    // virtual int code() const {return 1;}

    // const std::vector<double> &getMasses() const; 

    // void setInitial(int inid1, int inid2) 
    // {
    //   id1 = inid1; 
    //   id2 = inid2; 
    // }

    // int getDim() const {return dim;}

    // int getNIOParticles() const {return nexternal;}


    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 4; 
    static const int nprocesses = 1;
    static const int ncomb = 16; 

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
    int events_per_iter;

    // print verbose info
    // bool m_verbose; 

    // print debug info
    bool m_debug; 

    // static const int nwavefuncs = 6; 
    // static const int namplitudes = 2; 
    // static const int ncomb = 16; 
    // static const int wrows = 6; 
    // static const int nioparticles = 4;

    // Kokkos::complex<double> * * amp; 


    // Pointer to the model parameters
    Parameters_sm * pars; 

    // vector with external particle masses
    // std::vector<double> mME; 

    // Initial particle ids
    // int id1, id2; 

};
