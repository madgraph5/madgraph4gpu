//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <cmath> 
#include <thrust/complex.h> 


using namespace std; 

namespace MG5_sm 
{
__device__ void oxxxxx(double p[4], double fmass, int nhel, int nsf, 
thrust::complex<double> fo[6]); 

__device__ void sxxxxx(double p[4], int nss, thrust::complex<double> sc[3]); 

__device__ void ixxxxx(double p[4], double fmass, int nhel, int nsf, 
thrust::complex<double> fi[6]); 

__device__ void txxxxx(double p[4], double tmass, int nhel, int nst, 
thrust::complex<double> fi[18]); 

__device__ void vxxxxx(double p[4], double vmass, int nhel, int nsv, 
thrust::complex<double> v[6]); 

__device__ void FFV1_0(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    V3[], const __restrict__ thrust::complex<double> COUP, __restrict__
    thrust::complex<double> * vertex);

__device__ void FFV1P0_3(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    COUP, const __restrict__ double M3, const __restrict__ double W3,
    __restrict__ thrust::complex<double> V3[]);

__device__ void FFV2_0(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    V3[], const __restrict__ thrust::complex<double> COUP, __restrict__
    thrust::complex<double> * vertex);

__device__ void FFV2_3(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    COUP, const __restrict__ double M3, const __restrict__ double W3,
    __restrict__ thrust::complex<double> V3[]);

__device__ void FFV4_0(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    V3[], const __restrict__ thrust::complex<double> COUP, __restrict__
    thrust::complex<double> * vertex);

__device__ void FFV4_3(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    COUP, const __restrict__ double M3, const __restrict__ double W3,
    __restrict__ thrust::complex<double> V3[]);

__device__ void FFV2_4_0(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    V3[], const __restrict__ thrust::complex<double> COUP1, const __restrict__
    thrust::complex<double> COUP2, __restrict__ thrust::complex<double> *
    vertex);

__device__ void FFV2_4_3(thrust::complex<double> F1[], const __restrict__
    thrust::complex<double> F2[], const __restrict__ thrust::complex<double>
    COUP1, const __restrict__ thrust::complex<double> COUP2, const __restrict__
    double M3, const __restrict__ double W3, __restrict__
    thrust::complex<double> V3[]);

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
