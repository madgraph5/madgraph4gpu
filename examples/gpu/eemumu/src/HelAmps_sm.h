//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <cmath>

#include <thrust/complex.h>

namespace gMG5_sm {

__global__ void sigmaKin(cudaPitchedPtr tp, double *meDevPtr, size_t mePitch,
                         bool debug, bool verbose);

__device__ void matrix_1_epem_mupmum(double &matrix,
                                     thrust::complex<double> amp[2]);

__device__ void calculate_wavefunctions(int ihel, char *dps, size_t dpt,
                                        thrust::complex<double> amp[2],
                                        bool debug = false,
                                        bool verbose = false);

__device__ void oxxxxx(double p[4], double fmass, int nhel, int nsf,
                       thrust::complex<double> fo[6]);

__device__ void ixxxxx(double p[4], double fmass, int nhel, int nsf,
                       thrust::complex<double> fi[6]);

__device__ void FFV2_4_3(thrust::complex<double> F1[],
                         thrust::complex<double> F2[],
                         thrust::complex<double> COUP1,
                         thrust::complex<double> COUP2, double M3, double W3,
                         thrust::complex<double> V3[]);

__device__ void FFV1P0_3(thrust::complex<double> F1[],
                         thrust::complex<double> F2[],
                         thrust::complex<double> COUP, double M3, double W3,
                         thrust::complex<double> V3[]);

__device__ void
FFV2_4_0(thrust::complex<double> F1[], thrust::complex<double> F2[],
         thrust::complex<double> V3a[], thrust::complex<double> V3b[],
         thrust::complex<double> COUP0, thrust::complex<double> COUP1,
         thrust::complex<double> COUP2, thrust::complex<double> vertex[]);

} // namespace gMG5_sm

#endif // HelAmps_sm_H
