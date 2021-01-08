//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <CL/sycl.hpp>
#include <cmath>
#include <complex> 

namespace MG5_sm 
{
void oxxxxx(double p[4], double fmass, int nhel, int nsf, 
std::complex<double> fo[6]); 

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, 
std::complex<double> fi[6]); 

void txxxxx(double p[4], double tmass, int nhel, int nst, 
std::complex<double> fi[18]); 

void vxxxxx(double p[4], double vmass, int nhel, int nsv, 
std::complex<double> v[6]); 

void FFV1_0(std::complex<double> F1[], const
    std::complex<double> F2[], const std::complex<double> V3[], const
    std::complex<double> COUP, std::complex<double> * vertex);

void FFV1_1(std::complex<double> F2[], const
    std::complex<double> V3[], const std::complex<double> COUP, const
    double M1, const double W1, std::complex<double> F1[]);

void FFV1_2(std::complex<double> F1[], const
    std::complex<double> V3[], const std::complex<double> COUP, const
    double M2, const double W2, std::complex<double> F2[]);

void FFV1P0_3(std::complex<double> F1[], const
    std::complex<double> F2[], const std::complex<double> COUP, const
    double M3, const double W3, std::complex<double> V3[]);

void VVVV3_0(std::complex<double> V1[], const
    std::complex<double> V2[], const std::complex<double> V3[], const
    std::complex<double> V4[], const std::complex<double> COUP,
    std::complex<double> * vertex);

void VVVV3P0_1(std::complex<double> V2[], const
    std::complex<double> V3[], const std::complex<double> V4[], const
    std::complex<double> COUP, const double M1, const double W1,
    std::complex<double> V1[]);

void VVVV1_0(std::complex<double> V1[], const
    std::complex<double> V2[], const std::complex<double> V3[], const
    std::complex<double> V4[], const std::complex<double> COUP,
    std::complex<double> * vertex);

void VVVV1P0_1(std::complex<double> V2[], const
    std::complex<double> V3[], const std::complex<double> V4[], const
    std::complex<double> COUP, const double M1, const double W1,
    std::complex<double> V1[]);

void VVVV4_0(std::complex<double> V1[], const
    std::complex<double> V2[], const std::complex<double> V3[], const
    std::complex<double> V4[], const std::complex<double> COUP,
    std::complex<double> * vertex);

void VVVV4P0_1(std::complex<double> V2[], const
    std::complex<double> V3[], const std::complex<double> V4[], const
    std::complex<double> COUP, const double M1, const double W1,
    std::complex<double> V1[]);

void VVV1_0(std::complex<double> V1[], const
    std::complex<double> V2[], const std::complex<double> V3[], const
    std::complex<double> COUP, std::complex<double> * vertex);

void VVV1P0_1(std::complex<double> V2[], const
    std::complex<double> V3[], const std::complex<double> COUP, const
    double M1, const double W1, std::complex<double> V1[]);

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
