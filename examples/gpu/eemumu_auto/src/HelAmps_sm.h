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
__device__ void oxxxxx(const double p[3], double fmass, int nhel, int nsf, 
thrust::complex<double> fo[6]); 
__device__ void omzxxx(const double p[3], int nhel, int nsf, 
thrust::complex<double> fo[6]); 
__device__ void opzxxx(const double p[3], int nhel, int nsf, 
thrust::complex<double> fo[6]); 
__device__ void oxzxxx(const double p[3], int nhel, int nsf, 
thrust::complex<double> fo[6]); 

__device__ void sxxxxx(const double p[3], int nss, thrust::complex<double>
    sc[3]);

__device__ void ixxxxx(const double p[3], double fmass, int nhel, int nsf, 
thrust::complex<double> fi[6]); 
__device__ void imzxxx(const double p[3], int nhel, int nsf, 
thrust::complex<double> fi[6]); 
__device__ void ipzxxx(const double p[3], int nhel, int nsf, 
thrust::complex<double> fi[6]); 
__device__ void ixzxxx(const double p[3], int nhel, int nsf, 
thrust::complex<double> fi[6]); 

// __device__ void txxxxx(double p[4], double tmass, int nhel, int nst,
// thrust::complex<double> fi[18]);

__device__ void vxxxxx(const double p[3], double vmass, int nhel, int nsv, 
thrust::complex<double> v[6]); 

__device__ void FFV1_0(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> V3[], const
    thrust::complex<double> COUP, thrust::complex<double> * vertex);

__device__ void FFV1P0_3(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> COUP, const
    double M3, const double W3, thrust::complex<double> V3[]);

__device__ void FFV2_0(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> V3[], const
    thrust::complex<double> COUP, thrust::complex<double> * vertex);

__device__ void FFV2_3(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> COUP, const
    double M3, const double W3, thrust::complex<double> V3[]);

__device__ void FFV4_0(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> V3[], const
    thrust::complex<double> COUP, thrust::complex<double> * vertex);

__device__ void FFV4_3(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> COUP, const
    double M3, const double W3, thrust::complex<double> V3[]);

__device__ void FFV2_4_0(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> V3[], const
    thrust::complex<double> COUP1, const thrust::complex<double> COUP2,
    thrust::complex<double> * vertex);

__device__ void FFV2_4_3(const thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> COUP1, const
    thrust::complex<double> COUP2, const double M3, const double W3,
    thrust::complex<double> V3[]);

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
