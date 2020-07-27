//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include "mg5Complex.h"
#include <cmath>

using namespace std;

namespace MG5_sm
{
__device__ void oxxxxx(double p[4], double fmass, int nhel, int nsf,
mg5Complex fo[6]);

__device__ void sxxxxx(double p[4], int nss, mg5Complex sc[3]);

__device__ void ixxxxx(double p[4], double fmass, int nhel, int nsf,
mg5Complex fi[6]);

__device__ void txxxxx(double p[4], double tmass, int nhel, int nst,
mg5Complex fi[18]);

__device__ void vxxxxx(double p[4], double vmass, int nhel, int nsv,
mg5Complex v[6]);

__device__ void FFV1_0(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex V3[], const
    mg5Complex COUP, mg5Complex * vertex);

__device__ void FFV1P0_3(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex COUP, const
    double M3, const double W3, mg5Complex V3[]);

__device__ void FFV2_0(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex V3[], const
    mg5Complex COUP, mg5Complex * vertex);

__device__ void FFV2_3(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex COUP, const
    double M3, const double W3, mg5Complex V3[]);

__device__ void FFV4_0(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex V3[], const
    mg5Complex COUP, mg5Complex * vertex);

__device__ void FFV4_3(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex COUP, const
    double M3, const double W3, mg5Complex V3[]);

__device__ void FFV2_4_0(mg5Complex F1[], const
    mg5Complex F2[], const mg5Complex V3[], const
    mg5Complex COUP1, const mg5Complex COUP2,
    mg5Complex * vertex);

__device__ void FFV2_4_3(mg5Complex F1[], const mg5Complex F2[],
                         const mg5Complex COUP1, const mg5Complex COUP2,
                         const double M3, const double W3, mg5Complex V3[]);

}  // end namespace MG5_sm

#endif  // HelAmps_sm_H
