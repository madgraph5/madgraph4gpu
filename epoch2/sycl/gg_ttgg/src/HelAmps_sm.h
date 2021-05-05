//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#ifdef SYCL_LANGUAGE_VERSION
#include <CL/sycl.hpp>
#endif
#include <cmath>
#include "mgOnGpuConfig.h"
#include "mgOnGpuTypes.h"

using namespace std; 

namespace MG5_sm 
{
// #ifdef __CUDACC__
// 
// __device__ void oxxxxx(const fptype p[3], const fptype& fmass, const int&
// nhel, const int& nsf,
// cxtype fo[6]);
// __device__ void omzxxx(const fptype p[3],  const int& nhel, const int& nsf,
// cxtype fo[6]);
// __device__ void opzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fo[6]);
// __device__ void oxzxxx(const fptype p[3],  const int& nhel, const int& nsf,
// cxtype fo[6]);
// 
// __device__ void sxxxxx(const fptype p[3], const int& nss, cxtype sc[3]);
// 
// __device__ void ixxxxx(const fptype p[3], const fptype& fmass, const int&
// nhel, const int& nsf,
// cxtype fi[6]);
// __device__ void imzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fi[6]);
// __device__ void ipzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fi[6]);
// __device__ void ixzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fi[6]);
// 
// __device__ void txxxxx(fptype p[4], fptype tmass, int nhel, int nst,
// cxtype fi[18]);
// 
// __device__ void vxxxxx(const fptype p[3], const fptype& vmass, const int&
// nhel, const int& nsv,
// cxtype v[6]);
// 
// #else
// 
// __device__ void oxxxxx(const fptype p[3], const fptype& fmass, const int&
// nhel, const int& nsf, const int& ievt,
// cxtype fo[6]);
// __device__ void omzxxx(const fptype p[3],  const int& nhel, const int& nsf,
// cxtype fo[6]);
// __device__ void opzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fo[6]);
// __device__ void oxzxxx(const fptype p[3],  const int& nhel, const int& nsf,
// cxtype fo[6]);
// 
// __device__ void sxxxxx(const fptype p[3], const int& nss, cxtype sc[3]);
// 
// __device__ void ixxxxx(const fptype p[3], const fptype& fmass, const int&
// nhel, const int& nsf,
// cxtype fi[6]);
// __device__ void imzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fi[6]);
// __device__ void ipzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fi[6]);
// __device__ void ixzxxx(const fptype p[3], const int& nhel, const int& nsf,
// cxtype fi[6]);
// 
// __device__ void txxxxx(fptype p[4], fptype tmass, int nhel, int nst,
// cxtype fi[18]);
// 
// __device__ void vxxxxx(const fptype p[3], const fptype& vmass, const int&
// nhel, const int& nsv,
// cxtype v[6]);
// 
// 
// 
// #endif
// 
// 

void VVV1_0(const cxtype V1[], const cxtype V2[], const cxtype V3[],
    const cxtype COUP, cxtype * vertex);

void VVV1P0_1(const cxtype V2[], const cxtype V3[], const cxtype
    COUP, const fptype M1, const fptype W1, cxtype V1[]);

void VVVV1_0(const cxtype V1[], const cxtype V2[], const cxtype
    V3[], const cxtype V4[], const cxtype COUP, cxtype * vertex);

void VVVV1P0_1(const cxtype V2[], const cxtype V3[], const cxtype
    V4[], const cxtype COUP, const fptype M1, const fptype W1, cxtype V1[]);

void VVVV3_0(const cxtype V1[], const cxtype V2[], const cxtype
    V3[], const cxtype V4[], const cxtype COUP, cxtype * vertex);

void VVVV3P0_1(const cxtype V2[], const cxtype V3[], const cxtype
    V4[], const cxtype COUP, const fptype M1, const fptype W1, cxtype V1[]);

void FFV1_0(const cxtype F1[], const cxtype F2[], const cxtype V3[],
    const cxtype COUP, cxtype * vertex);

void FFV1_1(const cxtype F2[], const cxtype V3[], const cxtype COUP,
    const fptype M1, const fptype W1, cxtype F1[]);

void FFV1_2(const cxtype F1[], const cxtype V3[], const cxtype COUP,
    const fptype M2, const fptype W2, cxtype F2[]);

void FFV1P0_3(const cxtype F1[], const cxtype F2[], const cxtype
    COUP, const fptype M3, const fptype W3, cxtype V3[]);

void VVVV4_0(const cxtype V1[], const cxtype V2[], const cxtype
    V3[], const cxtype V4[], const cxtype COUP, cxtype * vertex);

void VVVV4P0_1(const cxtype V2[], const cxtype V3[], const cxtype
    V4[], const cxtype COUP, const fptype M1, const fptype W1, cxtype V1[]);

}  // end namespace MG5_sm

#endif // HelAmps_sm_H
