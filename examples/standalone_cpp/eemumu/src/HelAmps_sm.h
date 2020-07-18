//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <extras.h>

namespace MG5_sm {
double Sgn(double e, double f);

void oxxxxx(double p[4], double fmass, int nhel, int nsf,
            extras::complex fo[6]);

void sxxxxx(double p[4], int nss, extras::complex sc[3]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf,
            extras::complex fi[6]);

void txxxxx(double p[4], double tmass, int nhel, int nst,
            extras::complex fi[18]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv,
            extras::complex v[6]);

void FFV1_0(extras::complex F1[], extras::complex F2[],
            extras::complex V3[], extras::complex COUP,
            extras::complex &vertex);

void FFV2_3(extras::complex F1[], extras::complex F2[],
            extras::complex COUP, double M3, double W3,
            extras::complex V3[]);
void FFV2_4_3(extras::complex F1[], extras::complex F2[],
              extras::complex COUP1, extras::complex COUP2, double M3,
              double W3, extras::complex V3[]);

void FFV1P0_3(extras::complex F1[], extras::complex F2[],
              extras::complex COUP, double M3, double W3,
              extras::complex V3[]);

void FFV4_0(extras::complex F1[], extras::complex F2[],
            extras::complex V3[], extras::complex COUP,
            extras::complex &vertex);

void FFV4_3(extras::complex F1[], extras::complex F2[],
            extras::complex COUP, double M3, double W3,
            extras::complex V3[]);

void FFV2_0(extras::complex F1[], extras::complex F2[],
            extras::complex V3[], extras::complex COUP,
            extras::complex &vertex);
void FFV2_4_0(extras::complex F1[], extras::complex F2[],
              extras::complex V3[], extras::complex COUP1,
              extras::complex COUP2, extras::complex &vertex);

} // end namespace MG5_sm

#endif // HelAmps_sm_H
