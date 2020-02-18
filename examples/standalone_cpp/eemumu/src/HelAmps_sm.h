//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <cmath>
#include <complex>

using namespace std;

namespace MG5_sm {
double Sgn(double e, double f);

void oxxxxx(double p[4], double fmass, int nhel, int nsf,
            std::complex<double> fo[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf,
            std::complex<double> fi[6]);

void txxxxx(double p[4], double tmass, int nhel, int nst,
            std::complex<double> fi[18]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv,
            std::complex<double> v[6]);

void FFV1_0(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> V3[], std::complex<double> COUP,
            std::complex<double> &vertex);

void FFV2_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, double M3, double W3,
            std::complex<double> V3[]);
void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
              std::complex<double> COUP1, std::complex<double> COUP2, double M3,
              double W3, std::complex<double> V3[]);

void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
              std::complex<double> COUP, double M3, double W3,
              std::complex<double> V3[]);

void FFV4_0(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> V3[], std::complex<double> COUP,
            std::complex<double> &vertex);

void FFV4_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, double M3, double W3,
            std::complex<double> V3[]);

void FFV2_0(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> V3[], std::complex<double> COUP,
            std::complex<double> &vertex);
void FFV2_4_0(std::complex<double> F1[], std::complex<double> F2[],
              std::complex<double> V3[], std::complex<double> COUP1,
              std::complex<double> COUP2, std::complex<double> &vertex);

} // end namespace MG5_sm

#endif // HelAmps_sm_H
