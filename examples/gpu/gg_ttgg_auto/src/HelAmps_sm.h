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

__device__ void FFV1_0(thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> V3[], const
    thrust::complex<double> COUP, thrust::complex<double> * vertex);

__device__ void FFV1_1(thrust::complex<double> F2[], const
    thrust::complex<double> V3[], const thrust::complex<double> COUP, const
    double M1, const double W1, thrust::complex<double> F1[]);

__device__ void FFV1_2(thrust::complex<double> F1[], const
    thrust::complex<double> V3[], const thrust::complex<double> COUP, const
    double M2, const double W2, thrust::complex<double> F2[]);

__device__ void FFV1P0_3(thrust::complex<double> F1[], const
    thrust::complex<double> F2[], const thrust::complex<double> COUP, const
    double M3, const double W3, thrust::complex<double> V3[]);

__device__ void VVVV3_0(thrust::complex<double> V1[], const
    thrust::complex<double> V2[], const thrust::complex<double> V3[], const
    thrust::complex<double> V4[], const thrust::complex<double> COUP,
    thrust::complex<double> * vertex);

__device__ void VVVV3P0_1(thrust::complex<double> V2[], const
    thrust::complex<double> V3[], const thrust::complex<double> V4[], const
    thrust::complex<double> COUP, const double M1, const double W1,
    thrust::complex<double> V1[]);

__device__ void VVVV1_0(thrust::complex<double> V1[], const
    thrust::complex<double> V2[], const thrust::complex<double> V3[], const
    thrust::complex<double> V4[], const thrust::complex<double> COUP,
    thrust::complex<double> * vertex);

__device__ void VVVV1P0_1(thrust::complex<double> V2[], const
    thrust::complex<double> V3[], const thrust::complex<double> V4[], const
    thrust::complex<double> COUP, const double M1, const double W1,
    thrust::complex<double> V1[]);

__device__ void VVVV4_0(thrust::complex<double> V1[], const
    thrust::complex<double> V2[], const thrust::complex<double> V3[], const
    thrust::complex<double> V4[], const thrust::complex<double> COUP,
    thrust::complex<double> * vertex);

__device__ void VVVV4P0_1(thrust::complex<double> V2[], const
    thrust::complex<double> V3[], const thrust::complex<double> V4[], const
    thrust::complex<double> COUP, const double M1, const double W1,
    thrust::complex<double> V1[]);

__device__ void VVV1_0(thrust::complex<double> V1[], const
    thrust::complex<double> V2[], const thrust::complex<double> V3[], const
    thrust::complex<double> COUP, thrust::complex<double> * vertex);

__device__ void VVV1P0_1(thrust::complex<double> V2[], const
    thrust::complex<double> V3[], const thrust::complex<double> COUP, const
    double M1, const double W1, thrust::complex<double> V1[]);

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
