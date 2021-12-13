//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#include <string>
#include "Parameters_sm.h"
#include "Kokkos_Core.hpp"
// #define THRUST_COMPLEX 1
#ifdef THRUST_COMPLEX
  #include <thrust/complex.h>
  template<typename T>
  using complex = thrust::complex<T>;
#else
  #include "Kokkos_Complex.hpp"
  template<typename T>
  using complex = Kokkos::complex<T>;
#endif
#ifndef __CUDACC__
#include <cmath>
using std::min;
using std::max;
using std::sqrt;
using std::abs;
#endif


template<typename T>
KOKKOS_FUNCTION void ixxxxx(const T& pvec, const double fmass, 
  const int nhel, const int nsf, complex<double> fi[6])
{
  
  const double p0 = sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) + fmass * fmass);
  fi[0] = complex<double> (-p0 * nsf, -pvec(3) * nsf);
  fi[1] = complex<double> (-pvec(1) * nsf, -pvec(2) * nsf);
  if (fmass != 0.0)
  {
    double pp = min(p0, sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3)));
    if (pp == 0.0)
    {
      double sqm[2] = {sqrt(abs(fmass)), 0.};
      sqm[1] = (fmass < 0) ? -sqm[0] : sqm[0];
      const int ip = (1 + (nhel * nsf))/2;
      const int im = (1 - (nhel * nsf))/2;
      fi[2] = complex<double>(ip * sqm[ip], 0.);
      fi[3] = complex<double>(im * nsf * sqm[ip], 0.);
      fi[4] = complex<double>(ip * nsf * sqm[im], 0.);
      fi[5] = complex<double>(im * sqm[im], 0.);
    }
    else
    {
      const double sf[2] = { (1 + nsf + (1 - nsf) * (nhel * nsf)) * 0.5,
                             (1 + nsf - (1 - nsf) * (nhel * nsf)) * 0.5 };
      double omega[2] = { sqrt(p0 + pp), 0. };
      omega[1] = fmass / omega[0];
      const int ip = (1 + (nhel * nsf))/2;
      const int im = (1 - (nhel * nsf))/2;
      const double sfomega[2] = { sf[0] * omega[ip], sf[1] * omega[im] };
      const double pp3 = max(pp + pvec(3), 0.0);
      
      const complex<double> chi[2] = { complex<double> (sqrt(pp3 * 0.5/pp), 0),
                                               ( pp3 == 0. ?
                                                 complex<double> (-(nhel * nsf), 0) :
                                                 complex<double> ((nhel * nsf) * pvec(1), pvec(2)) 
                                                                  / sqrt(2.0 * pp * pp3) ) };
      fi[2] = sfomega[0] * chi[im];
      fi[3] = sfomega[0] * chi[ip];
      fi[4] = sfomega[1] * chi[im];
      fi[5] = sfomega[1] * chi[ip];
    }
  }
  else
  { 
    const double sqp0p3 = (pvec(1) == 0.0 and pvec(2) == 0.0 and pvec(3) < 0.0) ? 0. : sqrt(max(p0 + pvec(3), 0.0)) * nsf;

    const complex<double> chi[2] = {complex<double> (sqp0p3, 0.0),
          (sqp0p3 == 0.0) ? complex<double> (-nhel * sqrt(2.0 * p0), 0.0) :
                            complex<double> ((nhel * nsf) * pvec(1), pvec(2))/sqp0p3 };
    if ((nhel * nsf) == 1)
    {
      fi[2] = complex<double> (0.0, 0.0);
      fi[3] = complex<double> (0.0, 0.0);
      fi[4] = chi[0];
      fi[5] = chi[1];
    }
    else
    {
      fi[2] = chi[1];
      fi[3] = chi[0];
      fi[4] = complex<double> (0.0, 0.0);
      fi[5] = complex<double> (0.0, 0.0);
    }
  }
  return;
}


template<typename T>
KOKKOS_FUNCTION void vxxxxx(const T& pvec, const double vmass, 
  const int nhel, const int nsv, complex<double> vc[6]) 
{
  const double p0 = sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) + vmass * vmass);
  const double sqh = sqrt(0.5);
  const double hel = nhel;
  const int nsvahl = nsv * abs(hel);
  vc[0] = complex<double> ( p0 * nsv, pvec(3) * nsv );
  vc[1] = complex<double> ( pvec(1) * nsv, pvec(2) * nsv );
  if ( vmass != 0. )
  {
    const double pt2 = (pvec(1) * pvec(1)) + (pvec(2) * pvec(2));
    const double pp = min(p0, sqrt(pt2 + (pvec(3) * pvec(3))));
    const double pt = min(pp, sqrt(pt2));
    const double hel0 = 1. - abs( hel );
    if ( pp == 0. )
    {
      vc[2] = complex<double> ( 0., 0. );
      vc[3] = complex<double> ( -hel * sqh, 0. );
      vc[4] = complex<double> ( 0., nsvahl * sqh );
      vc[5] = complex<double> ( hel0, 0. );
    }
    else
    {
      const double emp = p0 / ( vmass * pp );
      vc[2] = complex<double> ( hel0 * pp / vmass, 0. );
      vc[5] = 
      complex<double> ( hel0 * pvec(3) * emp + hel * pt / pp * sqh, 0. );
      if ( pt != 0. )
      {
        const double pzpt = pvec(3) / ( pp * pt ) * sqh * hel;
        vc[3] = complex<double> ( hel0 * pvec(1) * emp - pvec(1) * pzpt, -nsvahl * pvec(2) / pt * sqh);
        vc[4] = complex<double> ( hel0 * pvec(2) * emp - pvec(2) * pzpt, nsvahl * pvec(1) / pt * sqh);
      }
      else
      {
        vc[3] = complex<double> ( -hel * sqh, 0. );
        vc[4] = complex<double> (0., nsvahl * (pvec(3) < 0  ? -sqh : sqh ) );
      }
    }
  }
  else
  {
    const double pt = sqrt((pvec(1) * pvec(1)) + (pvec(2) * pvec(2)));
    vc[2] = complex<double> (0., 0. );
    vc[5] = complex<double> (hel * pt / p0 * sqh, 0. );
    if ( pt != 0. )
    {
      const double pzpt = pvec(3) / (p0 * pt) * sqh * hel;
      vc[3] = complex<double> ( -pvec(1) * pzpt, -nsv * pvec(2) / pt * sqh);
      vc[4] = complex<double> ( -pvec(2) * pzpt, nsv * pvec(1) / pt * sqh);
    }
    else
    {
      vc[3] = complex<double> ( -hel * sqh, 0. );
      vc[4] = complex<double> ( 0., nsv * ( pvec(3) < 0 ? -sqh : sqh ) );
    }
  }
  return;
}


template<typename T>
KOKKOS_FUNCTION void oxxxxx(const T& pvec, const double fmass, 
  const int nhel, const int nsf, complex<double> fo[6]) 
{
  
  const double p0 = sqrt(pvec(1) * pvec(1) + pvec(2) * pvec(2) + pvec(3) * pvec(3) + fmass * fmass);

  fo[0] = complex<double> (p0 * nsf, pvec(3) * nsf);
  fo[1] = complex<double> (pvec(1) * nsf, pvec(2) * nsf);
  if (fmass != 0.)
  {
    const double pp = min(p0, sqrt((pvec(1) * pvec(1)) + (pvec(2) * pvec(2)) + (pvec(3) * pvec(3))));
    if (pp == 0.)
    {
      double sqm[2] = {sqrt(abs(fmass)), 0.};
      sqm[1] = fmass < 0 ? -sqm[0] : sqm[0];
      const int ip = -((1 - (nhel * nsf))/2) * nhel;
      const int im = (1 + (nhel * nsf))/2 * nhel;
      fo[2] = complex<double> (im * sqm[abs(ip)], 0.);
      fo[3] = complex<double> (ip * nsf * sqm[abs(ip)], 0.);
      fo[4] = complex<double> (im * nsf * sqm[abs(im)], 0.);
      fo[5] = complex<double> (ip * sqm[abs(im)], 0.);
    }
    else
    {
      const double sf[]= {double(1 + nsf + (1 - nsf) * (nhel * nsf)) * 0.5,
                          double(1 + nsf - (1 - nsf) * (nhel * nsf)) * 0.5};
      double omega[2] = {sqrt(p0 + pp), 0.};
      omega[1] = fmass/omega[0];
      const int ip = (1 + (nhel * nsf))/2;
      const int im = (1 - (nhel * nsf))/2;
      double sfomeg[2] = { sf[0] * omega[ip], sf[1] * omega[im]};
      const double pp3 = max(pp + pvec(3), 0.);
      const complex<double> chi[2] = { complex<double> (sqrt(pp3 * 0.5/pp), 0.),
                    pp3 == 0. ? complex<double> (-(nhel * nsf), 0.00) :
                                complex<double> ((nhel * nsf) * pvec(1), -pvec(2))/sqrt(2.0 * pp * pp3)
        };
      fo[2] = sfomeg[1] * chi[im];
      fo[3] = sfomeg[1] * chi[ip];
      fo[4] = sfomeg[0] * chi[im];
      fo[5] = sfomeg[0] * chi[ip];
    }
  }
  else
  {
    double sqp0p3 = ((pvec(1) == 0.00) and (pvec(2) == 0.00) and (pvec(3) < 0.00)) ?
              0. : sqrt(max(p0 + pvec(3), 0.00)) * nsf;

    const complex<double> chi[2] = { complex<double> (sqp0p3, 0.00),
            sqp0p3 == 0. ? complex<double> (-nhel, 0.) * sqrt(2. * p0) :
                           complex<double> ((nhel * nsf) * pvec(1), -pvec(2))/sqp0p3
    };

    if ((nhel * nsf) == 1)
    {
      fo[2] = chi[0];
      fo[3] = chi[1];
      fo[4] = complex<double> (0., 0.);
      fo[5] = complex<double> (0., 0.);
    }
    else
    {
      fo[2] = complex<double> (0., 0.);
      fo[3] = complex<double> (0., 0.);
      fo[4] = chi[1];
      fo[5] = chi[0];
    }
  }
  return;
}


KOKKOS_FUNCTION void VVVV3_0(const complex<double> V1[], 
    const complex<double> V2[], const complex<double> V3[], 
    const complex<double> V4[], const complex<double> COUP, 
    complex<double> * vertex)
{
  const complex<double> cI = complex<double> (0., 1.);
  const complex<double> TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  const complex<double> TMP9 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  const complex<double> TMP11 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  const complex<double> TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  (*vertex) = COUP * (-cI * (TMP6 * TMP9) + cI * (TMP1 * TMP11)); 
}


KOKKOS_FUNCTION void VVVV3P0_1(const complex<double> V2[], 
    const complex<double> V3[], const complex<double> V4[], 
    const complex<double> COUP, const double M1, const double W1, 
    complex<double> V1[])
{
  const complex<double> cI = complex<double> (0., 1.);
  V1[0] = +V2[0] + V3[0] + V4[0];
  V1[1] = +V2[1] + V3[1] + V4[1];
  const double P1[4] = {
    -V1[0].real(),
    -V1[1].real(),
    -V1[1].imag(),
    -V1[0].imag()};
  const complex<double> TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  const complex<double> TMP11 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  const complex<double> denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) 
                                              - (P1[3] * P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (TMP6 * V4[2]) + cI * (V2[2] * TMP11)); 
  V1[3] = denom * (-cI * (TMP6 * V4[3]) + cI * (V2[3] * TMP11)); 
  V1[4] = denom * (-cI * (TMP6 * V4[4]) + cI * (V2[4] * TMP11)); 
  V1[5] = denom * (-cI * (TMP6 * V4[5]) + cI * (V2[5] * TMP11));
}


KOKKOS_FUNCTION void VVVV1_0(const complex<double> V1[], 
    const complex<double> V2[], const complex<double> V3[], 
    const complex<double> V4[], const complex<double> COUP, 
    complex<double> * vertex)
{
  const complex<double> cI = complex<double> (0., 1.);
  const complex<double> TMP3 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  const complex<double> TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  const complex<double> TMP9 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  const complex<double> TMP10 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  (*vertex) = COUP * (-cI * (TMP6 * TMP9) + cI * (TMP3 * TMP10)); 
}


KOKKOS_FUNCTION void VVVV1P0_1(const complex<double> V2[], const
complex<double> V3[], const complex<double> V4[], const
complex<double> COUP, const double M1, const double W1, 
complex<double> V1[])
{
  const complex<double> cI = complex<double> (0., 1.);
  V1[0] = +V2[0] + V3[0] + V4[0];
  V1[1] = +V2[1] + V3[1] + V4[1];
  const double P1[] = {
    -V1[0].real(),
    -V1[1].real(), 
    -V1[1].imag(), 
    -V1[0].imag(),
  };
  const complex<double> TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  const complex<double> TMP10 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  const complex<double> denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) 
                                              - (P1[3] * P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (TMP6 * V4[2]) + cI * (V3[2] * TMP10)); 
  V1[3] = denom * (-cI * (TMP6 * V4[3]) + cI * (V3[3] * TMP10)); 
  V1[4] = denom * (-cI * (TMP6 * V4[4]) + cI * (V3[4] * TMP10)); 
  V1[5] = denom * (-cI * (TMP6 * V4[5]) + cI * (V3[5] * TMP10)); 
}


KOKKOS_FUNCTION void FFV1_0(const complex<double> F1[], const
complex<double> F2[], const complex<double> V3[], const
complex<double> COUP, complex<double> * vertex)
{
  const complex<double> cI = complex<double> (0., 1.);
 
  const complex<double> TMP12 = 
      (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] * (V3[2] + V3[5])))));
  (*vertex) = COUP * - cI * TMP12; 
}


KOKKOS_FUNCTION void FFV1_1(const complex<double> F2[], const
complex<double> V3[], const complex<double> COUP, const
double M1, const double W1, complex<double> F1[])
{
  const complex<double> cI = complex<double> (0., 1.);
  F1[0] = +F2[0] + V3[0];
  F1[1] = +F2[1] + V3[1];
  const double P1[] = {
    -F1[0].real(),
    -F1[1].real(), 
    -F1[1].imag(), 
    -F1[0].imag(),
  };
  constexpr double one( 1. );
  const complex<double> denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - 
                              (P1[2] * P1[2]) - (P1[3] * P1[3]) - M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (-V3[2] + V3[5]) + (P1[1] * (V3[3] -
      cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) + V3[4]) + P1[3] * (-V3[2] +
      V3[5])))) + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * (-one) *
      (V3[2] + V3[5]) + (P1[2] * (-one) * (+cI * (V3[2] + V3[5])) + P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] *
      (V3[3] + cI * (V3[4])))));
  F1[3] = denom * (-cI) * (F2[2] * (P1[0] * (-V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      (-one) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (-V3[3] + cI * (V3[4])) + F2[5] *
      (-V3[2] + V3[5]))));
  F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (-V3[3]
      + cI * (V3[4])) + (P1[2] * (-one) * (+cI * (V3[3]) + V3[4]) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (-V3[2] + V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (-V3[2] + V3[5]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (-V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * (-one) * (+cI * (V3[2] + V3[5])) + P1[3] *
      (-V3[3] + cI * (V3[4]))))) + (F2[5] * (P1[0] * (-V3[2] + V3[5]) + (P1[1]
      * (V3[3] + cI * (V3[4])) + (P1[2] * (-cI * (V3[3]) + V3[4]) + P1[3] *
      (-V3[2] + V3[5])))) + M1 * (F2[2] * (-V3[3] + cI * (V3[4])) + F2[3] *
      (V3[2] + V3[5]))));
}


KOKKOS_FUNCTION void FFV1_2(const complex<double> F1[], const
complex<double> V3[], const complex<double> COUP, const
double M2, const double W2, complex<double> F2[])
{
  const complex<double> cI = complex<double> (0., 1.);
  constexpr double one( 1. );
  F2[0] = +F1[0] + V3[0];
  F2[1] = +F1[1] + V3[1];
  const double P2[] = {
    -F2[0].real(),
    -F2[1].real(), 
    -F2[1].imag(), 
    -F2[0].imag(),
  };
  const complex<double> denom = COUP/((P2[0] * P2[0]) - (P2[1] * P2[1]) - 
                  (P2[2] * P2[2]) - (P2[3] * P2[3]) - M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (-one) *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (-V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2]) - cI * (V3[5])) + P2[3] *
      (-V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] *
      (-V3[3] + cI * (V3[4])))));
  F2[3] = denom * (-cI) * (F1[2] * (P2[0] * (-one) * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (-V3[2] + V3[5]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) + V3[4]) + P2[3] *
      (-V3[2] + V3[5])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] *
      (V3[2] + V3[5]))));
  F2[4] = denom * (-cI) * (F1[4] * (P2[0] * (-V3[2] + V3[5]) + (P2[1] * (V3[3]
      + cI * (V3[4])) + (P2[2] * (-cI * (V3[3]) + V3[4]) + P2[3] * (-V3[2] +
      V3[5])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * (-one) *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * (-one) * (V3[2] + V3[5]) + F1[3] * (-V3[3] +
      cI * (V3[4])))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * (-one) * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] - V3[5]) + (P2[2] * (+cI * (V3[2]) - cI * (V3[5])) +
      P2[3] * (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) +
      (P2[1] * (-V3[3] + cI * (V3[4])) + (P2[2] * (-one) * (+cI * (V3[3]) +
      V3[4]) - P2[3] * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI *
      (V3[4])) + F1[3] * (V3[2] - V3[5]))));
}


KOKKOS_FUNCTION void FFV1P0_3(complex<double> F1[], const
complex<double> F2[], const complex<double> COUP, const
double M3, const double W3, complex<double> V3[])
{
  const complex<double> cI = complex<double> (0., 1.);
 
  V3[0] = +F1[0] + F2[0];
  V3[1] = +F1[1] + F2[1];
  const double P3[] = {
    -V3[0].real(),
    -V3[1].real(), 
    -V3[1].imag(), 
    -V3[0].imag(),
  };
  const complex<double> denom = COUP/((P3[0] * P3[0]) - (P3[1] * P3[1]) - 
                (P3[2] * P3[2]) - (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]);
  V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] + F1[4] * F2[3] + F1[5] * F2[2]);
  V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] * F2[4] + F1[4] * F2[3]));
  V3[5] = denom * (-cI) * (-F1[2] * F2[4] - F1[5] * F2[3] + F1[3] * F2[5] + F1[4] * F2[2]);
}


KOKKOS_FUNCTION void VVVV4_0(const complex<double> V1[], const
complex<double> V2[], const complex<double> V3[], const
complex<double> V4[], const complex<double> COUP, 
complex<double> * vertex)
{
  const complex<double> cI = complex<double> (0., 1.);
 
  const complex<double> TMP3 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  const complex<double> TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  const complex<double> TMP10 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  const complex<double> TMP11 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  (*vertex) = COUP * (-cI * (TMP3 * TMP10) + cI * (TMP1 * TMP11)); 
}


KOKKOS_FUNCTION void VVVV4P0_1(complex<double> V2[], const
complex<double> V3[], const complex<double> V4[], const
complex<double> COUP, const double M1, const double W1, 
complex<double> V1[])
{
  const complex<double> cI = complex<double> (0., 1.);
  
  V1[0] = +V2[0] + V3[0] + V4[0];
  V1[1] = +V2[1] + V3[1] + V4[1];
  const double P1[] = {
    -V1[0].real(),
    -V1[1].real(), 
    -V1[1].imag(), 
    -V1[0].imag(),
  };
  const complex<double> TMP10 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  const complex<double> TMP11 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  const complex<double> denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - 
                                  (P1[2] * P1[2]) - (P1[3] * P1[3]) - M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V3[2] * TMP10) + cI * (V2[2] * TMP11)); 
  V1[3] = denom * (-cI * (V3[3] * TMP10) + cI * (V2[3] * TMP11)); 
  V1[4] = denom * (-cI * (V3[4] * TMP10) + cI * (V2[4] * TMP11)); 
  V1[5] = denom * (-cI * (V3[5] * TMP10) + cI * (V2[5] * TMP11)); 
}


KOKKOS_FUNCTION void VVV1_0(const complex<double> V1[], const
complex<double> V2[], const complex<double> V3[], const
complex<double> COUP, complex<double> * vertex)
{
  const complex<double> cI = complex<double> (0., 1.);
  const double P1[] = {
    V1[0].real(),
    V1[1].real(), 
    V1[1].imag(), 
    V1[0].imag(),
  };
  const double P2[] = {
    V2[0].real(),
    V2[1].real(), 
    V2[1].imag(), 
    V2[0].imag(),
  };
  const double P3[] = {
    V3[0].real(),
    V3[1].real(), 
    V3[1].imag(), 
    V3[0].imag(),
  };
  const complex<double> TMP3 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  const complex<double> TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  const complex<double> TMP4 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  const complex<double> TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  const complex<double> TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  const complex<double> TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  const complex<double> TMP8 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  const complex<double> TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  const complex<double> TMP5 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  (*vertex) = COUP * (TMP1 * (-cI * (TMP0) + cI * (TMP2)) + (TMP3 * (+cI *
      (TMP4) - cI * (TMP5)) + TMP6 * (-cI * (TMP7) + cI * (TMP8))));
}


KOKKOS_FUNCTION void VVV1P0_1(const complex<double> V2[],
    const complex<double> V3[], 
    const complex<double> COUP, 
    const double M1, const double W1, 
    complex<double> V1[])
{
  const complex<double> cI = complex<double> (0., 1.);
  
  const double P2[] = {
    V2[0].real(),
    V2[1].real(), 
    V2[1].imag(), 
    V2[0].imag(),
  };
  const double P3[] = {
    V3[0].real(),
    V3[1].real(), 
    V3[1].imag(), 
    V3[0].imag(),
  };
  
  V1[0] = +V2[0] + V3[0];
  V1[1] = +V2[1] + V3[1];
  
  const double P1[] = {
    -V1[0].real(),
    -V1[1].real(), 
    -V1[1].imag(), 
    -V1[0].imag(),
  };

  const complex<double> TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  const complex<double> TMP4 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  const complex<double> TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  const complex<double> TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  const complex<double> TMP5 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  const complex<double> denom = COUP/((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
      P1[3]) - M1 * (M1 - cI * W1));

  V1[2] = denom * (TMP6 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[2] * (+cI * (TMP4) - cI * (TMP5))));
  V1[3] = denom * (TMP6 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[3] * (+cI * (TMP4) - cI * (TMP5))));
  V1[4] = denom * (TMP6 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[4] * (+cI * (TMP4) - cI * (TMP5))));
  V1[5] = denom * (TMP6 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[5] * (+cI * (TMP4) - cI * (TMP5))));
}


template <typename hel_t, typename mom_t, typename ipd_t, typename ipc_t>
KOKKOS_FUNCTION void calculate_wavefunctions(
    const hel_t& cHel,
    const mom_t& local_mom,
    const ipd_t& cIPD,
    const ipc_t& cIPC,
    double& matrix)
{
  complex<double> amp;
  // Calculate wavefunctions for all processes
  complex<double> w[26][6];
  constexpr int ncolor = 24;
  complex<double> jamp[ncolor];
  const complex<double> const_complex(0,1);
#ifdef THRUST_COMPLEX
  for(int i = 0; i < ncolor; i++ )
  {
    jamp[i] = complex<double>(0., 0.); 
  }
#endif

  vxxxxx(Kokkos::subview(local_mom,0,Kokkos::ALL), 0., cHel(0), -1, w[0]);
  vxxxxx(Kokkos::subview(local_mom,1,Kokkos::ALL), 0., cHel(1), -1, w[1]);
  oxxxxx(Kokkos::subview(local_mom,2,Kokkos::ALL), cIPD(0), cHel(2), +1, w[2]);
  ixxxxx(Kokkos::subview(local_mom,3,Kokkos::ALL), cIPD(0), cHel(3), -1, w[3]);
  vxxxxx(Kokkos::subview(local_mom,4,Kokkos::ALL), 0., cHel(4), +1, w[4]);
  vxxxxx(Kokkos::subview(local_mom,5,Kokkos::ALL), 0., cHel(5), +1, w[5]);
  
  VVV1P0_1(w[0], w[1], cIPC(0), 0., 0., w[6]);
  FFV1P0_3(w[3], w[2], cIPC(1), 0., 0., w[7]);

  // Amplitude(s) for diagram number 1
  VVVV1_0(w[6], w[7], w[4], w[5], cIPC(2),&amp);
  auto a = +const_complex * amp;
  auto b = -a;
  jamp[0] += a;
  jamp[1] += b;
  jamp[6] += b;
  jamp[7] += a;
  jamp[16] += b;
  jamp[17] += a;
  jamp[22] += a;
  jamp[23] += b;
  VVVV3_0(w[6], w[7], w[4], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[6] += b;
  jamp[12] += b;
  jamp[14] += a;
  jamp[18] += b;
  jamp[20] += a;
  jamp[22] += a;
  jamp[23] += b;
  VVVV4_0(w[6], w[7], w[4], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[7] += b;
  jamp[12] += b;
  jamp[14] += a;
  jamp[16] += a;
  jamp[17] += b;
  jamp[18] += b;
  jamp[20] += a;
  VVV1P0_1(w[6], w[4], cIPC(0), 0., 0.,w[8]);

  // Amplitude(s) for diagram number 2
  VVV1_0(w[7], w[5], w[8], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[6] += b;
  jamp[12] += b;
  jamp[14] += a;
  jamp[18] += b;
  jamp[20] += a;
  jamp[22] += a;
  jamp[23] += b;
  VVV1P0_1(w[6], w[5], cIPC(0), 0., 0., w[9]);

  // Amplitude(s) for diagram number 3
  VVV1_0(w[7], w[4], w[9], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[7] += b;
  jamp[12] += b;
  jamp[14] += a;
  jamp[16] += a;
  jamp[17] += b;
  jamp[18] += b;
  jamp[20] += a;
  VVV1P0_1(w[4], w[5], cIPC(0), 0., 0., w[10]);

  // Amplitude(s) for diagram number 4
  VVV1_0(w[6], w[7], w[10], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[1] += b;
  jamp[6] += b;
  jamp[7] += a;
  jamp[16] += b;
  jamp[17] += a;
  jamp[22] += a;
  jamp[23] += b;
  FFV1_1(w[2], w[4], cIPC(1), cIPD(0), cIPD(1), w[11]);
  FFV1_2(w[3], w[6], cIPC(1), cIPD(0), cIPD(1), w[12]);

  // Amplitude(s) for diagram number 5
  FFV1_0(w[12], w[11], w[5], cIPC(1), &amp);
  jamp[16] += +const_complex * amp;
  jamp[17] += -const_complex * amp;

  // Amplitude(s) for diagram number 6
  FFV1_0(w[3], w[11], w[9], cIPC(1), &amp);
  jamp[12] += +amp;
  jamp[14] += -amp;
  jamp[16] += -amp;
  jamp[17] += +amp;
  FFV1_2(w[3], w[5], cIPC(1), cIPD(0), cIPD(1), w[13]);

  // Amplitude(s) for diagram number 7
  FFV1_0(w[13], w[11], w[6], cIPC(1), &amp);
  jamp[12] += +const_complex * amp;
  jamp[14] += -const_complex * amp;
  FFV1_1(w[2], w[5], cIPC(1), cIPD(0), cIPD(1), w[14]);

  // Amplitude(s) for diagram number 8
  FFV1_0(w[12], w[14], w[4], cIPC(1), &amp);
  jamp[22] += +const_complex * amp;
  jamp[23] += -const_complex * amp;

  // Amplitude(s) for diagram number 9
  FFV1_0(w[3], w[14], w[8], cIPC(1), &amp);
  jamp[18] += +amp;
  jamp[20] += -amp;
  jamp[22] += -amp;
  jamp[23] += +amp;
  FFV1_2(w[3], w[4], cIPC(1), cIPD(0), cIPD(1), w[15]);

  // Amplitude(s) for diagram number 10
  FFV1_0(w[15], w[14], w[6], cIPC(1), &amp);
  jamp[18] += +const_complex * amp;
  jamp[20] += -const_complex * amp;
  FFV1_1(w[2], w[6], cIPC(1), cIPD(0), cIPD(1), w[16]);

  // Amplitude(s) for diagram number 11
  FFV1_0(w[15], w[16], w[5], cIPC(1), &amp);
  jamp[1] += +const_complex * amp;
  jamp[7] += -const_complex * amp;

  // Amplitude(s) for diagram number 12
  FFV1_0(w[15], w[2], w[9], cIPC(1), &amp);
  jamp[1] += +amp;
  jamp[7] += -amp;
  jamp[18] += -amp;
  jamp[20] += +amp;

  // Amplitude(s) for diagram number 13
  FFV1_0(w[13], w[16], w[4], cIPC(1), &amp);
  jamp[0] += +const_complex * amp;
  jamp[6] += -const_complex * amp;

  // Amplitude(s) for diagram number 14
  FFV1_0(w[13], w[2], w[8], cIPC(1), &amp);
  jamp[0] += +amp;
  jamp[6] += -amp;
  jamp[12] += -amp;
  jamp[14] += +amp;

  // Amplitude(s) for diagram number 15
  FFV1_0(w[3], w[16], w[10], cIPC(1),&amp);
  jamp[0] += +amp;
  jamp[1] += -amp;
  jamp[6] += -amp;
  jamp[7] += +amp;

  // Amplitude(s) for diagram number 16
  FFV1_0(w[12], w[2], w[10], cIPC(1), &amp);
  jamp[16] += +amp;
  jamp[17] += -amp;
  jamp[22] += -amp;
  jamp[23] += +amp;
  FFV1_1(w[2], w[0], cIPC(1), cIPD(0), cIPD(1), w[12]);
  FFV1_2(w[3], w[1], cIPC(1), cIPD(0), cIPD(1), w[16]);
  FFV1_1(w[12], w[4], cIPC(1), cIPD(0), cIPD(1), w[8]);

  // Amplitude(s) for diagram number 17
  FFV1_0(w[16], w[8], w[5], cIPC(1), &amp);
  jamp[3] += -amp;
  FFV1_1(w[12], w[5], cIPC(1), cIPD(0), cIPD(1), w[9]);

  // Amplitude(s) for diagram number 18
  FFV1_0(w[16], w[9], w[4], cIPC(1), &amp);
  jamp[5] += -amp;

  // Amplitude(s) for diagram number 19
  FFV1_0(w[16], w[12], w[10], cIPC(1), &amp);
  jamp[3] += +const_complex * amp;
  jamp[5] += -const_complex * amp;
  VVV1P0_1(w[1], w[4], cIPC(0), 0., 0., w[6]);
  FFV1P0_3(w[3], w[12], cIPC(1), 0., 0., w[17]);

  // Amplitude(s) for diagram number 20
  VVV1_0(w[6], w[5], w[17], cIPC(0), &amp);
  jamp[0] += +amp;
  jamp[2] += -amp;
  jamp[4] += -amp;
  jamp[5] += +amp;

  // Amplitude(s) for diagram number 21
  FFV1_0(w[3], w[9], w[6], cIPC(1), &amp);
  jamp[4] += +const_complex * amp;
  jamp[5] += -const_complex * amp;

  // Amplitude(s) for diagram number 22
  FFV1_0(w[13], w[12], w[6], cIPC(1), &amp);
  jamp[0] += +const_complex * amp;
  jamp[2] += -const_complex * amp;
  VVV1P0_1(w[1], w[5], cIPC(0), 0., 0., w[18]);

  // Amplitude(s) for diagram number 23
  VVV1_0(w[18], w[4], w[17], cIPC(0), &amp);
  jamp[1] += +amp;
  jamp[2] += -amp;
  jamp[3] += +amp;
  jamp[4] += -amp;

  // Amplitude(s) for diagram number 24
  FFV1_0(w[3], w[8], w[18], cIPC(1), &amp);
  jamp[2] += +const_complex * amp;
  jamp[3] += -const_complex * amp;

  // Amplitude(s) for diagram number 25
  FFV1_0(w[15], w[12], w[18], cIPC(1), &amp);
  jamp[1] += +const_complex * amp;
  jamp[4] += -const_complex * amp;
  FFV1_1(w[12], w[1], cIPC(1), cIPD(0), cIPD(1), w[19]);

  // Amplitude(s) for diagram number 26
  FFV1_0(w[15], w[19], w[5], cIPC(1), &amp);
  jamp[1] += -amp;

  // Amplitude(s) for diagram number 27
  FFV1_0(w[15], w[9], w[1], cIPC(1), &amp);
  jamp[4] += -amp;

  // Amplitude(s) for diagram number 28
  FFV1_0(w[13], w[19], w[4], cIPC(1), &amp);
  jamp[0] += -amp;

  // Amplitude(s) for diagram number 29
  FFV1_0(w[13], w[8], w[1], cIPC(1), &amp);
  jamp[2] += -amp;

  // Amplitude(s) for diagram number 30
  FFV1_0(w[3], w[19], w[10], cIPC(1), &amp);
  jamp[0] += +const_complex * amp;
  jamp[1] += -const_complex * amp;

  // Amplitude(s) for diagram number 31
  VVV1_0(w[1], w[10], w[17], cIPC(0), &amp);
  jamp[0] += +amp;
  jamp[1] += -amp;
  jamp[3] += -amp;
  jamp[5] += +amp;
  VVVV1P0_1(w[1], w[4], w[5], cIPC(2), 0., 0., w[17]);
  VVVV3P0_1(w[1], w[4], w[5], cIPC(2), 0., 0., w[19]);
  VVVV4P0_1(w[1], w[4], w[5], cIPC(2), 0., 0., w[8]);

  // Amplitude(s) for diagram number 32
  FFV1_0(w[3], w[12], w[17], cIPC(1), &amp);
  jamp[0] += +amp;
  jamp[1] += -amp;
  jamp[3] += -amp;
  jamp[5] += +amp;
  FFV1_0(w[3], w[12], w[19], cIPC(1), &amp);
  jamp[1] += -amp;
  jamp[2] += +amp;
  jamp[3] += -amp;
  jamp[4] += +amp;
  FFV1_0(w[3], w[12], w[8], cIPC(1), &amp);
  jamp[0] += -amp;
  jamp[2] += +amp;
  jamp[4] += +amp;
  jamp[5] += -amp;
  FFV1_2(w[3], w[0], cIPC(1), cIPD(0), cIPD(1), w[12]);
  FFV1_1(w[2], w[1], cIPC(1), cIPD(0), cIPD(1), w[9]);
  FFV1_2(w[12], w[4], cIPC(1), cIPD(0), cIPD(1), w[20]);

  // Amplitude(s) for diagram number 33
  FFV1_0(w[20], w[9], w[5], cIPC(1), &amp);
  jamp[11] += -amp;
  FFV1_2(w[12], w[5], cIPC(1), cIPD(0), cIPD(1), w[21]);

  // Amplitude(s) for diagram number 34
  FFV1_0(w[21], w[9], w[4], cIPC(1), &amp);
  jamp[9] += -amp;

  // Amplitude(s) for diagram number 35
  FFV1_0(w[12], w[9], w[10], cIPC(1), &amp);
  jamp[9] += +const_complex * amp;
  jamp[11] += -const_complex * amp;
  FFV1P0_3(w[12], w[2], cIPC(1), 0., 0., w[22]);

  // Amplitude(s) for diagram number 36
  VVV1_0(w[6], w[5], w[22], cIPC(0), &amp);
  jamp[9] += +amp;
  jamp[15] += -amp;
  jamp[21] += -amp;
  jamp[23] += +amp;

  // Amplitude(s) for diagram number 37
  FFV1_0(w[21], w[2], w[6], cIPC(1), &amp);
  jamp[9] += +const_complex * amp;
  jamp[15] += -const_complex * amp;

  // Amplitude(s) for diagram number 38
  FFV1_0(w[12], w[14], w[6], cIPC(1), &amp);
  jamp[21] += +const_complex * amp;
  jamp[23] += -const_complex * amp;

  // Amplitude(s) for diagram number 39
  VVV1_0(w[18], w[4], w[22], cIPC(0), &amp);
  jamp[11] += +amp;
  jamp[15] += -amp;
  jamp[17] += +amp;
  jamp[21] += -amp;

  // Amplitude(s) for diagram number 40
  FFV1_0(w[20], w[2], w[18], cIPC(1), &amp);
  jamp[11] += +const_complex * amp;
  jamp[21] += -const_complex * amp;
  // Amplitude(s) for diagram number 41
  FFV1_0(w[12], w[11], w[18], cIPC(1), &amp);
  jamp[15] += +const_complex * amp;
  jamp[17] += -const_complex * amp;
  FFV1_2(w[12], w[1], cIPC(1), cIPD[0], cIPD[1], w[23]);
  // Amplitude(s) for diagram number 42
  FFV1_0(w[23], w[11], w[5], cIPC(1), &amp);
  jamp[17] += -amp;
  // Amplitude(s) for diagram number 43
  FFV1_0(w[21], w[11], w[1], cIPC(1), &amp);
  jamp[15] += -amp;
  // Amplitude(s) for diagram number 44
  FFV1_0(w[23], w[14], w[4], cIPC(1), &amp);
  jamp[23] += -amp;
  // Amplitude(s) for diagram number 45
  FFV1_0(w[20], w[14], w[1], cIPC(1), &amp);
  jamp[21] += -amp;
  // Amplitude(s) for diagram number 46
  FFV1_0(w[23], w[2], w[10], cIPC(1), &amp);
  jamp[17] += +const_complex * amp;
  jamp[23] += -const_complex * amp;
  // Amplitude(s) for diagram number 47
  VVV1_0(w[1], w[10], w[22], cIPC(0), &amp);
  jamp[9] += +amp;
  jamp[11] += -amp;
  jamp[17] += -amp;
  jamp[23] += +amp;
  // Amplitude(s) for diagram number 48
  FFV1_0(w[12], w[2], w[17], cIPC(1), &amp);
  jamp[9] += +amp;
  jamp[11] += -amp;
  jamp[17] += -amp;
  jamp[23] += +amp;
  FFV1_0(w[12], w[2], w[19], cIPC(1), &amp);
  jamp[11] += -amp;
  jamp[15] += +amp;
  jamp[17] += -amp;
  jamp[21] += +amp;
  FFV1_0(w[12], w[2], w[8], cIPC(1), &amp);
  jamp[9] += -amp;
  jamp[15] += +amp;
  jamp[21] += +amp;
  jamp[23] += -amp;
  VVV1P0_1(w[0], w[4], cIPC(0), 0., 0., w[12]);
  FFV1_2(w[3], w[12], cIPC(1), cIPD[0], cIPD[1], w[22]);
  // Amplitude(s) for diagram number 49
  FFV1_0(w[22], w[9], w[5], cIPC(1), &amp);
  jamp[10] += +const_complex * amp;
  jamp[11] += -const_complex * amp;
  VVV1P0_1(w[12], w[5], cIPC(0), 0., 0., w[23]);
  // Amplitude(s) for diagram number 50
  FFV1_0(w[3], w[9], w[23], cIPC(1), &amp);
  jamp[6] += +amp;
  jamp[8] += -amp;
  jamp[10] += -amp;
  jamp[11] += +amp;
  // Amplitude(s) for diagram number 51
  FFV1_0(w[13], w[9], w[12], cIPC(1), &amp);
  jamp[6] += +const_complex * amp;
  jamp[8] += -const_complex * amp;
  FFV1_1(w[2], w[12], cIPC(1), cIPD[0], cIPD[1], w[20]);
  // Amplitude(s) for diagram number 52
  FFV1_0(w[16], w[20], w[5], cIPC(1), &amp);
  jamp[3] += +const_complex * amp;
  jamp[13] += -const_complex * amp;
  // Amplitude(s) for diagram number 53
  FFV1_0(w[16], w[2], w[23], cIPC(1), &amp);
  jamp[3] += +amp;
  jamp[13] += -amp;
  jamp[19] += -amp;
  jamp[22] += +amp;
  // Amplitude(s) for diagram number 54
  FFV1_0(w[16], w[14], w[12], cIPC(1), &amp);
  jamp[19] += +const_complex * amp;
  jamp[22] += -const_complex * amp;
  // Amplitude(s) for diagram number 55
  FFV1_0(w[3], w[20], w[18], cIPC(1), &amp);
  jamp[2] += +amp;
  jamp[3] += -amp;
  jamp[12] += -amp;
  jamp[13] += +amp;
  // Amplitude(s) for diagram number 56
  FFV1_0(w[22], w[2], w[18], cIPC(1), &amp);
  jamp[10] += +amp;
  jamp[11] += -amp;
  jamp[20] += -amp;
  jamp[21] += +amp;
  // Amplitude(s) for diagram number 57
  VVV1_0(w[12], w[18], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[2] += b;
  jamp[3] += a;
  jamp[10] += a;
  jamp[11] += b;
  jamp[12] += a;
  jamp[13] += b;
  jamp[20] += b;
  jamp[21] += a;
  // Amplitude(s) for diagram number 58
  VVVV1_0(w[12], w[1], w[7], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[2] += a;
  jamp[6] += b;
  jamp[8] += a;
  jamp[12] += b;
  jamp[19] += b;
  jamp[20] += a;
  jamp[21] += b;
  jamp[22] += a;
  VVVV3_0(w[12], w[1], w[7], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[2] += a;
  jamp[3] += b;
  jamp[10] += b;
  jamp[11] += a;
  jamp[12] += b;
  jamp[13] += a;
  jamp[20] += a;
  jamp[21] += b;
  VVVV4_0(w[12], w[1], w[7], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[3] += b;
  jamp[6] += a;
  jamp[8] += b;
  jamp[10] += b;
  jamp[11] += a;
  jamp[13] += a;
  jamp[19] += a;
  jamp[22] += b;
  VVV1P0_1(w[12], w[1], cIPC(0), 0., 0., w[21]);
  // Amplitude(s) for diagram number 59
  VVV1_0(w[7], w[5], w[21], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[2] += a;
  jamp[6] += b;
  jamp[8] += a;
  jamp[12] += b;
  jamp[19] += b;
  jamp[20] += a;
  jamp[21] += b;
  jamp[22] += a;
  // Amplitude(s) for diagram number 60
  VVV1_0(w[1], w[7], w[23], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[3] += b;
  jamp[6] += a;
  jamp[8] += b;
  jamp[10] += b;
  jamp[11] += a;
  jamp[13] += a;
  jamp[19] += a;
  jamp[22] += b;
  // Amplitude(s) for diagram number 61
  FFV1_0(w[3], w[14], w[21], cIPC(1), &amp);
  jamp[19] += +amp;
  jamp[20] += -amp;
  jamp[21] += +amp;
  jamp[22] += -amp;
  // Amplitude(s) for diagram number 62
  FFV1_0(w[22], w[14], w[1], cIPC(1), &amp);
  jamp[20] += +const_complex * amp;
  jamp[21] += -const_complex * amp;
  // Amplitude(s) for diagram number 63
  FFV1_0(w[13], w[2], w[21], cIPC(1), &amp);
  jamp[2] += +amp;
  jamp[6] += -amp;
  jamp[8] += +amp;
  jamp[12] += -amp;
  // Amplitude(s) for diagram number 64
  FFV1_0(w[13], w[20], w[1], cIPC(1), &amp);
  jamp[2] += +const_complex * amp;
  jamp[12] += -const_complex * amp;
  VVV1P0_1(w[0], w[5], cIPC(0), 0., 0., w[20]);
  FFV1_2(w[3], w[20], cIPC(1), cIPD[0], cIPD[1], w[21]);
  // Amplitude(s) for diagram number 65
  FFV1_0(w[21], w[9], w[4], cIPC(1), &amp);
  jamp[8] += +const_complex * amp;
  jamp[9] += -const_complex * amp;
  VVV1P0_1(w[20], w[4], cIPC(0), 0., 0., w[22]);
  // Amplitude(s) for diagram number 66
  FFV1_0(w[3], w[9], w[22], cIPC(1), &amp);
  jamp[7] += +amp;
  jamp[8] += -amp;
  jamp[9] += +amp;
  jamp[10] += -amp;
  // Amplitude(s) for diagram number 67
  FFV1_0(w[15], w[9], w[20], cIPC(1), &amp);
  jamp[7] += +const_complex * amp;
  jamp[10] += -const_complex * amp;
  FFV1_1(w[2], w[20], cIPC(1), cIPD[0], cIPD[1], w[23]);
  // Amplitude(s) for diagram number 68
  FFV1_0(w[16], w[23], w[4], cIPC(1), &amp);
  jamp[5] += +const_complex * amp;
  jamp[19] += -const_complex * amp;
  // Amplitude(s) for diagram number 69
  FFV1_0(w[16], w[2], w[22], cIPC(1), &amp);
  jamp[5] += +amp;
  jamp[13] += -amp;
  jamp[16] += +amp;
  jamp[19] += -amp;
  // Amplitude(s) for diagram number 70
  FFV1_0(w[16], w[11], w[20], cIPC(1), &amp);
  jamp[13] += +const_complex * amp;
  jamp[16] += -const_complex * amp;
  // Amplitude(s) for diagram number 71
  FFV1_0(w[3], w[23], w[6], cIPC(1), &amp);
  jamp[4] += +amp;
  jamp[5] += -amp;
  jamp[18] += -amp;
  jamp[19] += +amp;
  // Amplitude(s) for diagram number 72
  FFV1_0(w[21], w[2], w[6], cIPC(1), &amp);
  jamp[8] += +amp;
  jamp[9] += -amp;
  jamp[14] += -amp;
  jamp[15] += +amp;
  // Amplitude(s) for diagram number 73
  VVV1_0(w[20], w[6], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[4] += b;
  jamp[5] += a;
  jamp[8] += a;
  jamp[9] += b;
  jamp[14] += b;
  jamp[15] += a;
  jamp[18] += a;
  jamp[19] += b;
  // Amplitude(s) for diagram number 74
  VVVV1_0(w[20], w[1], w[7], w[4], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[4] += a;
  jamp[7] += b;
  jamp[10] += a;
  jamp[13] += b;
  jamp[14] += a;
  jamp[15] += b;
  jamp[16] += a;
  jamp[18] += b;
  VVVV3_0(w[20], w[1], w[7], w[4], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[4] += a;
  jamp[5] += b;
  jamp[8] += b;
  jamp[9] += a;
  jamp[14] += a;
  jamp[15] += b;
  jamp[18] += b;
  jamp[19] += a;
  VVVV4_0(w[20], w[1], w[7], w[4], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[5] += b;
  jamp[7] += a;
  jamp[8] += b;
  jamp[9] += a;
  jamp[10] += b;
  jamp[13] += a;
  jamp[16] += b;
  jamp[19] += a;
  VVV1P0_1(w[20], w[1], cIPC(0), 0., 0., w[12]);
  // Amplitude(s) for diagram number 75
  VVV1_0(w[7], w[4], w[12], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[4] += a;
  jamp[7] += b;
  jamp[10] += a;
  jamp[13] += b;
  jamp[14] += a;
  jamp[15] += b;
  jamp[16] += a;
  jamp[18] += b;
  // Amplitude(s) for diagram number 76
  VVV1_0(w[1], w[7], w[22], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[5] += b;
  jamp[7] += a;
  jamp[8] += b;
  jamp[9] += a;
  jamp[10] += b;
  jamp[13] += a;
  jamp[16] += b;
  jamp[19] += a;
  // Amplitude(s) for diagram number 77
  FFV1_0(w[3], w[11], w[12], cIPC(1), &amp);
  jamp[13] += +amp;
  jamp[14] += -amp;
  jamp[15] += +amp;
  jamp[16] += -amp;
  // Amplitude(s) for diagram number 78
  FFV1_0(w[21], w[11], w[1], cIPC(1), &amp);
  jamp[14] += +const_complex * amp;
  jamp[15] += -const_complex * amp;
  // Amplitude(s) for diagram number 79
  FFV1_0(w[15], w[2], w[12], cIPC(1), &amp);
  jamp[4] += +amp;
  jamp[7] += -amp;
  jamp[10] += +amp;
  jamp[18] += -amp;
  // Amplitude(s) for diagram number 80
  FFV1_0(w[15], w[23], w[1], cIPC(1), &amp);
  jamp[4] += +const_complex * amp;
  jamp[18] += -const_complex * amp;
  FFV1_1(w[9], w[0], cIPC(1), cIPD[0], cIPD[1], w[23]);
  // Amplitude(s) for diagram number 81
  FFV1_0(w[15], w[23], w[5], cIPC(1), &amp);
  jamp[7] += -amp;
  FFV1_2(w[15], w[0], cIPC(1), cIPD[0], cIPD[1], w[12]);
  // Amplitude(s) for diagram number 82
  FFV1_0(w[12], w[9], w[5], cIPC(1), &amp);
  jamp[10] += -amp;
  // Amplitude(s) for diagram number 83
  FFV1_0(w[13], w[23], w[4], cIPC(1), &amp);
  jamp[6] += -amp;
  FFV1_2(w[13], w[0], cIPC(1), cIPD[0], cIPD[1], w[21]);
  // Amplitude(s) for diagram number 84
  FFV1_0(w[21], w[9], w[4], cIPC(1), &amp);
  jamp[8] += -amp;
  // Amplitude(s) for diagram number 85
  FFV1_0(w[3], w[23], w[10], cIPC(1), &amp);
  jamp[6] += +const_complex * amp;
  jamp[7] += -const_complex * amp;
  VVV1P0_1(w[0], w[10], cIPC(0), 0., 0., w[23]);
  // Amplitude(s) for diagram number 86
  FFV1_0(w[3], w[9], w[23], cIPC(1), &amp);
  jamp[6] += +amp;
  jamp[7] += -amp;
  jamp[9] += -amp;
  jamp[11] += +amp;
  FFV1_2(w[16], w[0], cIPC(1), cIPD[0], cIPD[1], w[22]);
  // Amplitude(s) for diagram number 87
  FFV1_0(w[22], w[11], w[5], cIPC(1), &amp);
  jamp[16] += -amp;
  FFV1_1(w[11], w[0], cIPC(1), cIPD[0], cIPD[1], w[20]);
  // Amplitude(s) for diagram number 88
  FFV1_0(w[16], w[20], w[5], cIPC(1), &amp);
  jamp[13] += -amp;
  // Amplitude(s) for diagram number 89
  FFV1_0(w[22], w[14], w[4], cIPC(1), &amp);
  jamp[22] += -amp;
  FFV1_1(w[14], w[0], cIPC(1), cIPD[0], cIPD[1], w[24]);
  // Amplitude(s) for diagram number 90
  FFV1_0(w[16], w[24], w[4], cIPC(1), &amp);
  jamp[19] += -amp;
  // Amplitude(s) for diagram number 91
  FFV1_0(w[22], w[2], w[10], cIPC(1), &amp);
  jamp[16] += +const_complex * amp;
  jamp[22] += -const_complex * amp;
  // Amplitude(s) for diagram number 92
  FFV1_0(w[16], w[2], w[23], cIPC(1), &amp);
  jamp[3] += +amp;
  jamp[5] += -amp;
  jamp[16] += -amp;
  jamp[22] += +amp;
  // Amplitude(s) for diagram number 93
  VVVV1_0(w[0], w[6], w[7], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[2] += b;
  jamp[8] += b;
  jamp[14] += a;
  jamp[18] += b;
  jamp[19] += a;
  jamp[21] += a;
  jamp[23] += b;
  VVVV3_0(w[0], w[6], w[7], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[2] += b;
  jamp[4] += b;
  jamp[5] += a;
  jamp[9] += b;
  jamp[15] += a;
  jamp[21] += a;
  jamp[23] += b;
  VVVV4_0(w[0], w[6], w[7], w[5], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[4] += b;
  jamp[5] += a;
  jamp[8] += a;
  jamp[9] += b;
  jamp[14] += b;
  jamp[15] += a;
  jamp[18] += a;
  jamp[19] += b;
  VVV1P0_1(w[0], w[6], cIPC(0), 0., 0., w[22]);
  // Amplitude(s) for diagram number 94
  VVV1_0(w[7], w[5], w[22], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[2] += b;
  jamp[8] += b;
  jamp[14] += a;
  jamp[18] += b;
  jamp[19] += a;
  jamp[21] += a;
  jamp[23] += b;
  VVV1P0_1(w[0], w[7], cIPC(0), 0., 0., w[25]);
  // Amplitude(s) for diagram number 95
  VVV1_0(w[6], w[5], w[25], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[2] += b;
  jamp[4] += b;
  jamp[5] += a;
  jamp[9] += b;
  jamp[15] += a;
  jamp[21] += a;
  jamp[23] += b;
  // Amplitude(s) for diagram number 96
  FFV1_0(w[3], w[14], w[22], cIPC(1), &amp);
  jamp[18] += +amp;
  jamp[19] += -amp;
  jamp[21] += -amp;
  jamp[23] += +amp;
  // Amplitude(s) for diagram number 97
  FFV1_0(w[3], w[24], w[6], cIPC(1), &amp);
  jamp[18] += +const_complex * amp;
  jamp[19] += -const_complex * amp;
  // Amplitude(s) for diagram number 98
  FFV1_0(w[13], w[2], w[22], cIPC(1), &amp);
  jamp[0] += +amp;
  jamp[2] += -amp;
  jamp[8] += -amp;
  jamp[14] += +amp;
  // Amplitude(s) for diagram number 99
  FFV1_0(w[21], w[2], w[6], cIPC(1), &amp);
  jamp[8] += +const_complex * amp;
  jamp[14] += -const_complex * amp;
  // Amplitude(s) for diagram number 100
  VVVV1_0(w[0], w[18], w[7], w[4], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[4] += b;
  jamp[10] += b;
  jamp[12] += b;
  jamp[13] += a;
  jamp[15] += a;
  jamp[17] += b;
  jamp[20] += a;
  VVVV3_0(w[0], w[18], w[7], w[4], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[2] += b;
  jamp[3] += a;
  jamp[4] += b;
  jamp[11] += b;
  jamp[15] += a;
  jamp[17] += b;
  jamp[21] += a;
  VVVV4_0(w[0], w[18], w[7], w[4], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[2] += b;
  jamp[3] += a;
  jamp[10] += a;
  jamp[11] += b;
  jamp[12] += a;
  jamp[13] += b;
  jamp[20] += b;
  jamp[21] += a;
  VVV1P0_1(w[0], w[18], cIPC(0), 0., 0., w[6]);
  // Amplitude(s) for diagram number 101
  VVV1_0(w[7], w[4], w[6], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[4] += b;
  jamp[10] += b;
  jamp[12] += b;
  jamp[13] += a;
  jamp[15] += a;
  jamp[17] += b;
  jamp[20] += a;
  // Amplitude(s) for diagram number 102
  VVV1_0(w[18], w[4], w[25], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[2] += b;
  jamp[3] += a;
  jamp[4] += b;
  jamp[11] += b;
  jamp[15] += a;
  jamp[17] += b;
  jamp[21] += a;
  // Amplitude(s) for diagram number 103
  FFV1_0(w[3], w[11], w[6], cIPC(1), &amp);
  jamp[12] += +amp;
  jamp[13] += -amp;
  jamp[15] += -amp;
  jamp[17] += +amp;
  // Amplitude(s) for diagram number 104
  FFV1_0(w[3], w[20], w[18], cIPC(1), &amp);
  jamp[12] += +const_complex * amp;
  jamp[13] += -const_complex * amp;
  // Amplitude(s) for diagram number 105
  FFV1_0(w[15], w[2], w[6], cIPC(1), &amp);
  jamp[1] += +amp;
  jamp[4] += -amp;
  jamp[10] += -amp;
  jamp[20] += +amp;
  // Amplitude(s) for diagram number 106
  FFV1_0(w[12], w[2], w[18], cIPC(1), &amp);
  jamp[10] += +const_complex * amp;
  jamp[20] += -const_complex * amp;
  // Amplitude(s) for diagram number 107
  VVVV1_0(w[0], w[1], w[7], w[10], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[1] += b;
  jamp[6] += b;
  jamp[7] += a;
  jamp[16] += b;
  jamp[17] += a;
  jamp[22] += a;
  jamp[23] += b;
  VVVV3_0(w[0], w[1], w[7], w[10], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[1] += b;
  jamp[3] += b;
  jamp[5] += a;
  jamp[9] += b;
  jamp[11] += a;
  jamp[17] += a;
  jamp[23] += b;
  VVVV4_0(w[0], w[1], w[7], w[10], cIPC(2), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[3] += b;
  jamp[5] += a;
  jamp[6] += a;
  jamp[7] += b;
  jamp[9] += b;
  jamp[11] += a;
  jamp[16] += a;
  jamp[22] += b;
  // Amplitude(s) for diagram number 108
  VVV1_0(w[1], w[10], w[25], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[1] += b;
  jamp[3] += b;
  jamp[5] += a;
  jamp[9] += b;
  jamp[11] += a;
  jamp[17] += a;
  jamp[23] += b;
  // Amplitude(s) for diagram number 109
  VVV1_0(w[1], w[7], w[23], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[3] += b;
  jamp[5] += a;
  jamp[6] += a;
  jamp[7] += b;
  jamp[9] += b;
  jamp[11] += a;
  jamp[16] += a;
  jamp[22] += b;
  // Amplitude(s) for diagram number 110
  FFV1_0(w[13], w[20], w[1], cIPC(1), &amp);
  jamp[12] += -amp;
  // Amplitude(s) for diagram number 111
  FFV1_0(w[21], w[11], w[1], cIPC(1), &amp);
  jamp[14] += -amp;
  // Amplitude(s) for diagram number 112
  FFV1_0(w[15], w[24], w[1], cIPC(1), &amp);
  jamp[18] += -amp;
  // Amplitude(s) for diagram number 113
  FFV1_0(w[12], w[14], w[1], cIPC(1), &amp);
  jamp[20] += -amp;
  VVVV1P0_1(w[0], w[1], w[4], cIPC(2), 0., 0., w[12]);
  VVVV3P0_1(w[0], w[1], w[4], cIPC(2), 0., 0., w[24]);
  VVVV4P0_1(w[0], w[1], w[4], cIPC(2), 0., 0., w[21]);
  // Amplitude(s) for diagram number 114
  VVV1_0(w[12], w[7], w[5], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[2] += b;
  jamp[8] += b;
  jamp[14] += a;
  jamp[18] += b;
  jamp[19] += a;
  jamp[21] += a;
  jamp[23] += b;
  VVV1_0(w[24], w[7], w[5], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[2] += b;
  jamp[6] += a;
  jamp[8] += b;
  jamp[12] += a;
  jamp[19] += a;
  jamp[20] += b;
  jamp[21] += a;
  jamp[22] += b;
  VVV1_0(w[21], w[7], w[5], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += b;
  jamp[6] += a;
  jamp[12] += a;
  jamp[14] += b;
  jamp[18] += a;
  jamp[20] += b;
  jamp[22] += b;
  jamp[23] += a;
  // Amplitude(s) for diagram number 115
  FFV1_0(w[3], w[14], w[12], cIPC(1), &amp);
  jamp[18] += +amp;
  jamp[19] += -amp;
  jamp[21] += -amp;
  jamp[23] += +amp;
  FFV1_0(w[3], w[14], w[24], cIPC(1), &amp);
  jamp[19] += -amp;
  jamp[20] += +amp;
  jamp[21] += -amp;
  jamp[22] += +amp;
  FFV1_0(w[3], w[14], w[21], cIPC(1), &amp);
  jamp[18] += -amp;
  jamp[20] += +amp;
  jamp[22] += +amp;
  jamp[23] += -amp;
  // Amplitude(s) for diagram number 116
  FFV1_0(w[13], w[2], w[12], cIPC(1), &amp);
  jamp[0] += +amp;
  jamp[2] += -amp;
  jamp[8] += -amp;
  jamp[14] += +amp;
  FFV1_0(w[13], w[2], w[24], cIPC(1), &amp);
  jamp[2] += -amp;
  jamp[6] += +amp;
  jamp[8] += -amp;
  jamp[12] += +amp;
  FFV1_0(w[13], w[2], w[21], cIPC(1), &amp);
  jamp[0] += -amp;
  jamp[6] += +amp;
  jamp[12] += +amp;
  jamp[14] += -amp;
  VVVV1P0_1(w[0], w[1], w[5], cIPC(2), 0., 0., w[21]);
  VVVV3P0_1(w[0], w[1], w[5], cIPC(2), 0., 0., w[13]);
  VVVV4P0_1(w[0], w[1], w[5], cIPC(2), 0., 0., w[24]);
  // Amplitude(s) for diagram number 117
  VVV1_0(w[21], w[7], w[4], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[4] += b;
  jamp[10] += b;
  jamp[12] += b;
  jamp[13] += a;
  jamp[15] += a;
  jamp[17] += b;
  jamp[20] += a;
  VVV1_0(w[13], w[7], w[4], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[4] += b;
  jamp[7] += a;
  jamp[10] += b;
  jamp[13] += a;
  jamp[14] += b;
  jamp[15] += a;
  jamp[16] += b;
  jamp[18] += a;
  VVV1_0(w[24], w[7], w[4], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += b;
  jamp[7] += a;
  jamp[12] += a;
  jamp[14] += b;
  jamp[16] += b;
  jamp[17] += a;
  jamp[18] += a;
  jamp[20] += b;
  // Amplitude(s) for diagram number 118
  FFV1_0(w[3], w[11], w[21], cIPC(1), &amp);
  jamp[12] += +amp;
  jamp[13] += -amp;
  jamp[15] += -amp;
  jamp[17] += +amp;
  FFV1_0(w[3], w[11], w[13], cIPC(1), &amp);
  jamp[13] += -amp;
  jamp[14] += +amp;
  jamp[15] += -amp;
  jamp[16] += +amp;
  FFV1_0(w[3], w[11], w[24], cIPC(1), &amp);
  jamp[12] += -amp;
  jamp[14] += +amp;
  jamp[16] += +amp;
  jamp[17] += -amp;
  // Amplitude(s) for diagram number 119
  FFV1_0(w[15], w[2], w[21], cIPC(1), &amp);
  jamp[1] += +amp;
  jamp[4] += -amp;
  jamp[10] += -amp;
  jamp[20] += +amp;
  FFV1_0(w[15], w[2], w[13], cIPC(1), &amp);
  jamp[4] += -amp;
  jamp[7] += +amp;
  jamp[10] += -amp;
  jamp[18] += +amp;
  FFV1_0(w[15], w[2], w[24], cIPC(1), &amp);
  jamp[1] += -amp;
  jamp[7] += +amp;
  jamp[18] += +amp;
  jamp[20] += -amp;
  VVVV1P0_1(w[0], w[4], w[5], cIPC(2), 0., 0., w[24]);
  VVVV3P0_1(w[0], w[4], w[5], cIPC(2), 0., 0., w[15]);
  VVVV4P0_1(w[0], w[4], w[5], cIPC(2), 0., 0., w[13]);
  // Amplitude(s) for diagram number 120
  FFV1_0(w[3], w[9], w[24], cIPC(1), &amp);
  jamp[6] += +amp;
  jamp[7] += -amp;
  jamp[9] += -amp;
  jamp[11] += +amp;
  FFV1_0(w[3], w[9], w[15], cIPC(1), &amp);
  jamp[7] += -amp;
  jamp[8] += +amp;
  jamp[9] += -amp;
  jamp[10] += +amp;
  FFV1_0(w[3], w[9], w[13], cIPC(1), &amp);
  jamp[6] += -amp;
  jamp[8] += +amp;
  jamp[10] += +amp;
  jamp[11] += -amp;
  // Amplitude(s) for diagram number 121
  FFV1_0(w[16], w[2], w[24], cIPC(1), &amp);
  jamp[3] += +amp;
  jamp[5] += -amp;
  jamp[16] += -amp;
  jamp[22] += +amp;
  FFV1_0(w[16], w[2], w[15], cIPC(1), &amp);
  jamp[5] += -amp;
  jamp[13] += +amp;
  jamp[16] += -amp;
  jamp[19] += +amp;
  FFV1_0(w[16], w[2], w[13], cIPC(1), &amp);
  jamp[3] += -amp;
  jamp[13] += +amp;
  jamp[19] += +amp;
  jamp[22] += -amp;
  // Amplitude(s) for diagram number 122
  VVV1_0(w[24], w[1], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[3] += b;
  jamp[5] += a;
  jamp[6] += a;
  jamp[7] += b;
  jamp[9] += b;
  jamp[11] += a;
  jamp[16] += a;
  jamp[22] += b;
  VVV1_0(w[15], w[1], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[5] += a;
  jamp[7] += b;
  jamp[8] += a;
  jamp[9] += b;
  jamp[10] += a;
  jamp[13] += b;
  jamp[16] += a;
  jamp[19] += b;
  VVV1_0(w[13], w[1], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[3] += a;
  jamp[6] += b;
  jamp[8] += a;
  jamp[10] += a;
  jamp[11] += b;
  jamp[13] += b;
  jamp[19] += b;
  jamp[22] += a;
  // Amplitude(s) for diagram number 123
  VVV1_0(w[0], w[17], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += b;
  jamp[1] += a;
  jamp[3] += a;
  jamp[5] += b;
  jamp[9] += a;
  jamp[11] += b;
  jamp[17] += b;
  jamp[23] += a;
  VVV1_0(w[0], w[19], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[1] += a;
  jamp[2] += b;
  jamp[3] += a;
  jamp[4] += b;
  jamp[11] += b;
  jamp[15] += a;
  jamp[17] += b;
  jamp[21] += a;
  VVV1_0(w[0], w[8], w[7], cIPC(0), &amp);
  a = +const_complex * amp;
  b = -a;
  jamp[0] += a;
  jamp[2] += b;
  jamp[4] += b;
  jamp[5] += a;
  jamp[9] += b;
  jamp[15] += a;
  jamp[21] += a;
  jamp[23] += b;
  
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
  for(int icol = 0;icol < ncolor; ++icol )
  {
    complex<double> ztemp;
    for(int jcol = 0;jcol < ncolor; ++jcol ){
      ztemp = ztemp + cf[icol][jcol] * jamp[jcol];
    }
    matrix = matrix + (ztemp * conj(jamp[icol])).real()/denom[icol];
  }
} // end calculate_wavefunction()

template <typename mom_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin_setup(
    const mom_t& momenta,
    const hel_t& cHel,
    const ipd_t& cIPD,
    const ipc_t& cIPC,
    igh_t& iGoodHel,
    ngh_t& nGoodHel,
    const int& ncomb,
    const int& league_size,
    const int& team_size) 
{
  Kokkos::View<int*,Kokkos::DefaultExecutionSpace> isGoodHel("isGoodHel",ncomb); // has to be constant index, but should be `ncomb`

  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(member_type team_member){
    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

    const int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
    double meHelSum[nprocesses] = {0};  // all zeros
    double meHelSumLast = 0;

    auto evt_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ihel = 0;ihel < ncomb;++ihel)
    {
      auto local_cHel = Kokkos::subview(cHel,ihel,Kokkos::ALL);
      calculate_wavefunctions(local_cHel, evt_mom, cIPD, cIPC, meHelSum[0]);
      if (meHelSum[0] != meHelSumLast)
      {
        isGoodHel(ihel) = true;
        meHelSumLast = meHelSum[0];
      }
    }
  });

  Kokkos::parallel_for(__func__,Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0,1),
  KOKKOS_LAMBDA(const int& i){
    for(int ihel=0;ihel < ncomb;++ihel){

      if(isGoodHel(ihel)){
        iGoodHel(nGoodHel(0)) = ihel;
        nGoodHel(0)++;
      }

    }
  });
}


//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
template <typename mom_t, typename out_t, typename hel_t, typename ipd_t, 
          typename ipc_t, typename igh_t, typename ngh_t>
void sigmaKin(const mom_t& momenta, out_t& allMEs, const hel_t& cHel,
    const ipd_t& cIPD, const ipc_t& cIPC, const igh_t& iGoodHel,
    const ngh_t& nGoodHel, const int& ncomb, const int& league_size,
    const int& team_size)
{
  using member_type = typename Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy( league_size, team_size );
  Kokkos::parallel_for(__func__,policy, 
  KOKKOS_LAMBDA(member_type team_member){

    const int ievt = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    
    // Denominators: spins, colors and identical particles
    constexpr int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
    constexpr int denominators[1] = {512};

    // Reset the "matrix elements" - running sums of |M|^2 over helicities for
    // the given event
    double meHelSum[nprocesses] = {0};  // all zeros

    // using precomputed good helicities
    auto local_mom = Kokkos::subview(momenta,ievt,Kokkos::ALL,Kokkos::ALL);
    for (int ighel = 0;ighel < nGoodHel(0);++ighel)
    {
      auto local_cHel = Kokkos::subview(cHel,iGoodHel(ighel),Kokkos::ALL);
      calculate_wavefunctions(local_cHel, local_mom, cIPD, cIPC, meHelSum[0]);
    }
    // Get the final |M|^2 as an average over helicities/colors of the running
    // sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf
    // 1d7c/Handout_4_2016-UZH.pdf]
    for (int iproc = 0;iproc < nprocesses; ++iproc)
    {
      meHelSum[iproc] /= denominators[iproc];
    }

    // Set the final average |M|^2 for this event in the output array for all
    // events
    for (int iproc = 0;iproc < nprocesses; ++iproc)
    {
      allMEs[iproc * nprocesses + ievt] = meHelSum[iproc];
    }
  });// end parallel for

}

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ g g WEIGHTED<=4 @1
//--------------------------------------------------------------------------

template <class ExecSpace>
class CPPProcess
{
  public:

    CPPProcess(int numiterations, int gpublocks, int gputhreads): 
        m_numiterations(numiterations), gpu_nblocks(gpublocks), 
        gpu_nthreads(gputhreads), 
        dim(gpu_nblocks * gpu_nthreads),
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
        {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1},{1, 1,
      1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
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

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 6;
    static const int nprocesses = 1;
    static const int ncomb = 64;

    Kokkos::View<int**,ExecSpace> cHel;
    typename Kokkos::View<int**,ExecSpace>::HostMirror hHel;

    Kokkos::View<double*,ExecSpace> cmME;
    typename Kokkos::View<double*,ExecSpace>::HostMirror hmME;

    Kokkos::View<complex<double>*,ExecSpace> cIPC;
    typename Kokkos::View<complex<double>*,ExecSpace>::HostMirror hIPC;

    Kokkos::View<double*,ExecSpace> cIPD;
    typename Kokkos::View<double*,ExecSpace>::HostMirror hIPD;

  private:
    int m_numiterations;
    // gpu variables
    int gpu_nblocks;
    int gpu_nthreads;
    int dim;  // gpu_nblocks * gpu_nthreads;
    Parameters_sm * pars;


};

// Constants for array limits
template <typename ExecSpace> const int CPPProcess<ExecSpace>::ninitial;
template <typename ExecSpace> const int CPPProcess<ExecSpace>::nprocesses;
template <typename ExecSpace> const int CPPProcess<ExecSpace>::ncomb;
template <typename ExecSpace> const int CPPProcess<ExecSpace>::nexternal;
