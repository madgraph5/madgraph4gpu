//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

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
#endif


namespace MG5_sm 
{
   template<typename T>
KOKKOS_FUNCTION void ixxxxx(const T& pvec, const double fmass, 
  const int nhel, const int nsf, complex<double> fi[6]);

template<typename T>
KOKKOS_FUNCTION void vxxxxx(const T& pvec, const double vmass, 
  const int nhel, const int nsv, complex<double> vc[6]);

template<typename T>
KOKKOS_FUNCTION void oxxxxx(const T& pvec, const double fmass, 
  const int nhel, const int nsf, complex<double> fo[6]);

KOKKOS_FUNCTION void VVVV3_0(const complex<double> V1[], 
    const complex<double> V2[], const complex<double> V3[], 
    const complex<double> V4[], const complex<double> COUP, 
    complex<double> * vertex);

KOKKOS_FUNCTION void VVVV3P0_1(const complex<double> V2[], 
    const complex<double> V3[], const complex<double> V4[], 
    const complex<double> COUP, const double M1, const double W1, 
    complex<double> V1[]);


KOKKOS_FUNCTION void VVVV1_0(const complex<double> V1[], 
    const complex<double> V2[], const complex<double> V3[], 
    const complex<double> V4[], const complex<double> COUP, 
    complex<double> * vertex);

KOKKOS_FUNCTION void VVVV1P0_1(const complex<double> V2[], const
complex<double> V3[], const complex<double> V4[], const
complex<double> COUP, const double M1, const double W1, 
complex<double> V1[]);

KOKKOS_FUNCTION void FFV1_0(const complex<double> F1[], const
complex<double> F2[], const complex<double> V3[], const
complex<double> COUP, complex<double> * vertex);

KOKKOS_FUNCTION void FFV1_1(const complex<double> F2[], const
complex<double> V3[], const complex<double> COUP, const
double M1, const double W1, complex<double> F1[]);

KOKKOS_FUNCTION void FFV1_2(const complex<double> F1[], const
complex<double> V3[], const complex<double> COUP, const
double M2, const double W2, complex<double> F2[]);

KOKKOS_FUNCTION void FFV1P0_3(complex<double> F1[], const
complex<double> F2[], const complex<double> COUP, const
double M3, const double W3, complex<double> V3[]);

KOKKOS_FUNCTION void VVVV4_0(const complex<double> V1[], const
complex<double> V2[], const complex<double> V3[], const
complex<double> V4[], const complex<double> COUP, 
complex<double> * vertex);

KOKKOS_FUNCTION void VVVV4P0_1(complex<double> V2[], const
complex<double> V3[], const complex<double> V4[], const
complex<double> COUP, const double M1, const double W1, 
complex<double> V1[]);

KOKKOS_FUNCTION void VVV1_0(const complex<double> V1[], const
complex<double> V2[], const complex<double> V3[], const
complex<double> COUP, complex<double> * vertex);

KOKKOS_FUNCTION void VVV1P0_1(const complex<double> V2[],
    const complex<double> V3[], 
    const complex<double> COUP, 
    const double M1, const double W1, 
    complex<double> V1[]);


}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
