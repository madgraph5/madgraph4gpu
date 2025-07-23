//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 3.6.3, 2025-06-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// GPU Template
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

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

struct ALOHAOBJ
{
    fptype p[4];
    cxtype W[4];
    int flv_index = 1;

    public:
        inline ALOHAOBJ() {};
}

__device__ void VVV1P0_1(ALOHAOBJ V2, ALOHAOBJ V3[], const cxtype
    COUP, const fptype M1, const fptype W1, ALOHAOBJ & V1);

__device__ void FFV1_0(ALOHAOBJ F1, ALOHAOBJ F2, ALOHAOBJ V3,
    const cxtype COUP, cxtype * vertex);

__device__ void FFV1_1(ALOHAOBJ F2, ALOHAOBJ V3, const cxtype COUP,
    const fptype M1, const fptype W1, ALOHAOBJ & F1);

__device__ void FFV1_2(ALOHAOBJ F1, ALOHAOBJ V3, const cxtype COUP,
    const fptype M2, const fptype W2, ALOHAOBJ & F2[]);

}  // end namespace MG5_sm

#endif // HelAmps_sm_H
