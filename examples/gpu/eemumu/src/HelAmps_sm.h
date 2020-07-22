//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <math.h>
#include <complex>

namespace gMG5_sm {

//void sigmaKin(cudaPitchedPtr tp, double *meDevPtr, size_t mePitch,
                         //bool debug, bool verbose);
void sigmaKin(const double[4][4], bool debug, bool verbose);

void matrix_1_epem_mupmum(double &matrix,
                                     std::complex<double> amp[2]);

void calculate_wavefunctions(int ihel, std::complex<double> amp[2],
                                        bool debug = false,
                                        bool verbose = false);

void oxxxxx(double p[4], double fmass, int nhel, int nsf,
                       std::complex<double> fo[6]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf,
                       std::complex<double> fi[6]);

void FFV2_4_3(std::complex<double> F1[],
                         std::complex<double> F2[],
                         std::complex<double> COUP1,
                         std::complex<double> COUP2, double M3, double W3,
                         std::complex<double> V3[]);

void FFV1P0_3(std::complex<double> F1[],
                         std::complex<double> F2[],
                         std::complex<double> COUP, double M3, double W3,
                         std::complex<double> V3[]);

void
FFV2_4_0(std::complex<double> F1[], std::complex<double> F2[],
         std::complex<double> V3a[], std::complex<double> V3b[],
         std::complex<double> COUP0, std::complex<double> COUP1,
         std::complex<double> COUP2, std::complex<double> vertex[]);

} // namespace gMG5_sm

#endif // HelAmps_sm_H
