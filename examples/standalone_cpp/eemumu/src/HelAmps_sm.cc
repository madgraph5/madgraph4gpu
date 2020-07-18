//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "HelAmps_sm.h"
#include <cmath>
#include <extras.h>
#include <cstdlib>
#include <iostream>

namespace MG5_sm {

void ixxxxx(double p[4], double fmass, int nhel, int nsf,
            extras::complex fi[6]) {
  extras::complex chi[2];
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2];
  int ip, im, nh;
  fi[0] = extras::complex(-p[0] * nsf, -p[3] * nsf);
  fi[1] = extras::complex(-p[1] * nsf, -p[2] * nsf);
  nh = nhel * nsf;
  if (fmass != 0.0) {
    pp = std::min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]));
    if (pp == 0.0) {
      sqm[0] = sqrt(std::abs(fmass));
      sqm[1] = Sgn(sqm[0], fmass);
      ip = (1 + nh) / 2;
      im = (1 - nh) / 2;
      fi[2] = ip * sqm[ip];
      fi[3] = im * nsf * sqm[ip];
      fi[4] = ip * nsf * sqm[im];
      fi[5] = im * sqm[im];
    } else {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5;
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5;
      omega[0] = sqrt(p[0] + pp);
      omega[1] = fmass / omega[0];
      ip = (1 + nh) / 2;
      im = (1 - nh) / 2;
      sfomega[0] = sf[0] * omega[ip];
      sfomega[1] = sf[1] * omega[im];
      pp3 = std::max(pp + p[3], 0.0);
      chi[0] = extras::complex(sqrt(pp3 * 0.5 / pp), 0);
      if (pp3 == 0.0) {
        chi[1] = extras::complex(-nh, 0);
      } else {
        chi[1] = extras::complex(nh * p[1], p[2]) / sqrt(2.0 * pp * pp3);
      }
      fi[2] = sfomega[0] * chi[im];
      fi[3] = sfomega[0] * chi[ip];
      fi[4] = sfomega[1] * chi[im];
      fi[5] = sfomega[1] * chi[ip];
    }
  } else {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0) {
      sqp0p3 = 0.0;
    } else {
      sqp0p3 = sqrt(std::max(p[0] + p[3], 0.0)) * nsf;
    }
    chi[0] = extras::complex(sqp0p3, 0.0);
    if (sqp0p3 == 0.0) {
      chi[1] = extras::complex(-nhel * sqrt(2.0 * p[0]), 0.0);
    } else {
      chi[1] = extras::complex(nh * p[1], p[2]) / sqp0p3;
    }
    if (nh == 1) {
      fi[2] = extras::complex(0.0, 0.0);
      fi[3] = extras::complex(0.0, 0.0);
      fi[4] = chi[0];
      fi[5] = chi[1];
    } else {
      fi[2] = chi[1];
      fi[3] = chi[0];
      fi[4] = extras::complex(0.0, 0.0);
      fi[5] = extras::complex(0.0, 0.0);
    }
  }
  return;
}

double Sgn(double a, double b) { return (b < 0) ? -std::abs(a) : std::abs(a); }

void txxxxx(double p[4], double tmass, int nhel, int nst,
            extras::complex tc[18]) {
  extras::complex ft[6][4], ep[4], em[4], e0[4];
  double pt, pt2, pp, pzpt, emp, sqh, sqs;
  int i, j;

  sqh = sqrt(0.5);
  sqs = sqrt(0.5 / 3);

  pt2 = p[1] * p[1] + p[2] * p[2];
  pp = std::min(p[0], sqrt(pt2 + p[3] * p[3]));
  pt = std::min(pp, sqrt(pt2));

  ft[4][0] = extras::complex(p[0] * nst, p[3] * nst);
  ft[5][0] = extras::complex(p[1] * nst, p[2] * nst);

  // construct eps+
  if (nhel >= 0) {
    if (pp == 0) {
      ep[0] = extras::complex(0, 0);
      ep[1] = extras::complex(-sqh, 0);
      ep[2] = extras::complex(0, nst * sqh);
      ep[3] = extras::complex(0, 0);
    } else {
      ep[0] = extras::complex(0, 0);
      ep[3] = extras::complex(pt / pp * sqh, 0);

      if (pt != 0) {
        pzpt = p[3] / (pp * pt) * sqh;
        ep[1] = extras::complex(-p[1] * pzpt, -nst * p[2] / pt * sqh);
        ep[2] = extras::complex(-p[2] * pzpt, nst * p[1] / pt * sqh);
      } else {
        ep[1] = extras::complex(-sqh, 0);
        ep[2] = extras::complex(0, nst * Sgn(sqh, p[3]));
      }
    }
  }

  // construct eps-
  if (nhel <= 0) {
    if (pp == 0) {
      em[0] = extras::complex(0, 0);
      em[1] = extras::complex(sqh, 0);
      em[2] = extras::complex(0, nst * sqh);
      em[3] = extras::complex(0, 0);
    } else {
      em[0] = extras::complex(0, 0);
      em[3] = extras::complex(-pt / pp * sqh, 0);

      if (pt != 0) {
        pzpt = -p[3] / (pp * pt) * sqh;
        em[1] = extras::complex(-p[1] * pzpt, -nst * p[2] / pt * sqh);
        em[2] = extras::complex(-p[2] * pzpt, nst * p[1] / pt * sqh);
      } else {
        em[1] = extras::complex(sqh, 0);
        em[2] = extras::complex(0, nst * Sgn(sqh, p[3]));
      }
    }
  }

  // construct eps0
  if (std::labs(nhel) <= 1) {
    if (pp == 0) {
      e0[0] = extras::complex(0, 0);
      e0[1] = extras::complex(0, 0);
      e0[2] = extras::complex(0, 0);
      e0[3] = extras::complex(1, 0);
    } else {
      emp = p[0] / (tmass * pp);
      e0[0] = extras::complex(pp / tmass, 0);
      e0[3] = extras::complex(p[3] * emp, 0);

      if (pt != 0) {
        e0[1] = extras::complex(p[1] * emp, 0);
        e0[2] = extras::complex(p[2] * emp, 0);
      } else {
        e0[1] = extras::complex(0, 0);
        e0[2] = extras::complex(0, 0);
      }
    }
  }

  if (nhel == 2) {
    for (j = 0; j < 4; j++) {
      for (i = 0; i < 4; i++)
        ft[i][j] = ep[i] * ep[j];
    }
  } else if (nhel == -2) {
    for (j = 0; j < 4; j++) {
      for (i = 0; i < 4; i++)
        ft[i][j] = em[i] * em[j];
    }
  } else if (tmass == 0) {
    for (j = 0; j < 4; j++) {
      for (i = 0; i < 4; i++)
        ft[i][j] = 0;
    }
  } else if (tmass != 0) {
    if (nhel == 1) {
      for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++)
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]);
      }
    } else if (nhel == 0) {
      for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++)
          ft[i][j] =
              sqs * (ep[i] * em[j] + em[i] * ep[j] + 2.0 * e0[i] * e0[j]);
      }
    } else if (nhel == -1) {
      for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++)
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]);
      }
    } else {
      std::cerr << "Invalid helicity in txxxxx.\n";
      std::exit(1);
    }
  }

  tc[0] = ft[4][0];
  tc[1] = ft[5][0];

  for (j = 0; j < 4; j++) {
    for (i = 0; i < 4; i++)
      tc[j * 4 + i + 2] = ft[j][i];
  }
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv,
            extras::complex vc[6]) {
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh;
  int nsvahl;
  sqh = sqrt(0.5);
  hel = double(nhel);
  nsvahl = nsv * std::abs(hel);
  pt2 = (p[1] * p[1]) + (p[2] * p[2]);
  pp = std::min(p[0], sqrt(pt2 + (p[3] * p[3])));
  pt = std::min(pp, sqrt(pt2));
  vc[0] = extras::complex(p[0] * nsv, p[3] * nsv);
  vc[1] = extras::complex(p[1] * nsv, p[2] * nsv);
  if (vmass != 0.0) {
    hel0 = 1.0 - std::abs(hel);
    if (pp == 0.0) {
      vc[2] = extras::complex(0.0, 0.0);
      vc[3] = extras::complex(-hel * sqh, 0.0);
      vc[4] = extras::complex(0.0, nsvahl * sqh);
      vc[5] = extras::complex(hel0, 0.0);
    } else {
      emp = p[0] / (vmass * pp);
      vc[2] = extras::complex(hel0 * pp / vmass, 0.0);
      vc[5] = extras::complex(hel0 * p[3] * emp + hel * pt / pp * sqh, 0.0);
      if (pt != 0.0) {
        pzpt = p[3] / (pp * pt) * sqh * hel;
        vc[3] = extras::complex(hel0 * p[1] * emp - p[1] * pzpt,
                                -nsvahl * p[2] / pt * sqh);
        vc[4] = extras::complex(hel0 * p[2] * emp - p[2] * pzpt,
                                nsvahl * p[1] / pt * sqh);
      } else {
        vc[3] = extras::complex(-hel * sqh, 0.0);
        vc[4] = extras::complex(0.0, nsvahl * Sgn(sqh, p[3]));
      }
    }
  } else {
    pp = p[0];
    pt = sqrt((p[1] * p[1]) + (p[2] * p[2]));
    vc[2] = extras::complex(0.0, 0.0);
    vc[5] = extras::complex(hel * pt / pp * sqh, 0.0);
    if (pt != 0.0) {
      pzpt = p[3] / (pp * pt) * sqh * hel;
      vc[3] = extras::complex(-p[1] * pzpt, -nsv * p[2] / pt * sqh);
      vc[4] = extras::complex(-p[2] * pzpt, nsv * p[1] / pt * sqh);
    } else {
      vc[3] = extras::complex(-hel * sqh, 0.0);
      vc[4] = extras::complex(0.0, nsv * Sgn(sqh, p[3]));
    }
  }
  return;
}

void sxxxxx(double p[4], int nss, extras::complex sc[3]) {
  sc[2] = extras::complex(1.00, 0.00);
  sc[0] = extras::complex(p[0] * nss, p[3] * nss);
  sc[1] = extras::complex(p[1] * nss, p[2] * nss);
  return;
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf,
            extras::complex fo[6]) {
  extras::complex chi[2];
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2];
  int nh, ip, im;
  fo[0] = extras::complex(p[0] * nsf, p[3] * nsf);
  fo[1] = extras::complex(p[1] * nsf, p[2] * nsf);
  nh = nhel * nsf;
  if (fmass != 0.000) {
    pp = std::min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3])));
    if (pp == 0.000) {
      sqm[0] = sqrt(std::abs(fmass));
      sqm[1] = Sgn(sqm[0], fmass);
      ip = -((1 - nh) / 2) * nhel;
      im = (1 + nh) / 2 * nhel;
      fo[2] = im * sqm[std::abs(ip)];
      fo[3] = ip * nsf * sqm[std::abs(ip)];
      fo[4] = im * nsf * sqm[std::abs(im)];
      fo[5] = ip * sqm[std::abs(im)];
    } else {
      pp = std::min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3])));
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5;
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5;
      omega[0] = sqrt(p[0] + pp);
      omega[1] = fmass / omega[0];
      ip = (1 + nh) / 2;
      im = (1 - nh) / 2;
      sfomeg[0] = sf[0] * omega[ip];
      sfomeg[1] = sf[1] * omega[im];
      pp3 = std::max(pp + p[3], 0.00);
      chi[0] = extras::complex(sqrt(pp3 * 0.5 / pp), 0.00);
      if (pp3 == 0.00) {
        chi[1] = extras::complex(-nh, 0.00);
      } else {
        chi[1] = extras::complex(nh * p[1], -p[2]) / sqrt(2.0 * pp * pp3);
      }
      fo[2] = sfomeg[1] * chi[im];
      fo[3] = sfomeg[1] * chi[ip];
      fo[4] = sfomeg[0] * chi[im];
      fo[5] = sfomeg[0] * chi[ip];
    }
  } else {
    if ((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00)) {
      sqp0p3 = 0.00;
    } else {
      sqp0p3 = sqrt(std::max(p[0] + p[3], 0.00)) * nsf;
    }
    chi[0] = extras::complex(sqp0p3, 0.00);
    if (sqp0p3 == 0.000) {
      chi[1] = extras::complex(-nhel, 0.00) * sqrt(2.0 * p[0]);
    } else {
      chi[1] = extras::complex(nh * p[1], -p[2]) / sqp0p3;
    }
    if (nh == 1) {
      fo[2] = chi[0];
      fo[3] = chi[1];
      fo[4] = extras::complex(0.00, 0.00);
      fo[5] = extras::complex(0.00, 0.00);
    } else {
      fo[2] = extras::complex(0.00, 0.00);
      fo[3] = extras::complex(0.00, 0.00);
      fo[4] = chi[1];
      fo[5] = chi[0];
    }
  }
  return;
}

void FFV1_0(extras::complex F1[], extras::complex F2[],
            extras::complex V3[], extras::complex COUP,
            extras::complex &vertex) {
  static extras::complex cI = extras::complex(0., 1.);
  extras::complex TMP2;
  TMP2 =
      (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
       (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
        (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
         F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * -cI * TMP2;
}

void FFV2_3(extras::complex F1[], extras::complex F2[],
            extras::complex COUP, double M3, double W3,
            extras::complex V3[]) {
  static extras::complex cI = extras::complex(0., 1.);
  extras::complex denom;
  extras::complex TMP1;
  double P3[4];
  double OM3;
  OM3 = 0.;
  if (M3 != 0.)
    OM3 = 1. / (M3 * M3);
  V3[0] = +F1[0] + F2[0];
  V3[1] = +F1[1] + F2[1];
  P3[0] = -V3[0].real();
  P3[1] = -V3[1].real();
  P3[2] = -V3[1].imag();
  P3[3] = -V3[0].imag();
  TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
          F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP / ((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) -
                  (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP1);
  V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 * TMP1);
  V3[4] = denom * (-cI) *
          (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2] * OM3 * TMP1);
  V3[5] = denom * (-cI) * (F1[3] * F2[5] - F1[2] * F2[4] - P3[3] * OM3 * TMP1);
}

void FFV2_4_3(extras::complex F1[], extras::complex F2[],
              extras::complex COUP1, extras::complex COUP2, double M3,
              double W3, extras::complex V3[]) {
  int i;
  extras::complex Vtmp[6];
  FFV2_3(F1, F2, COUP1, M3, W3, V3);
  FFV4_3(F1, F2, COUP2, M3, W3, Vtmp);
  i = 2;
  while (i < 6) {
    V3[i] = V3[i] + Vtmp[i];
    i++;
  }
}

void FFV1P0_3(extras::complex F1[], extras::complex F2[],
              extras::complex COUP, double M3, double W3,
              extras::complex V3[]) {
  static extras::complex cI = extras::complex(0., 1.);
  double P3[4];
  extras::complex denom;
  V3[0] = +F1[0] + F2[0];
  V3[1] = +F1[1] + F2[1];
  P3[0] = -V3[0].real();
  P3[1] = -V3[1].real();
  P3[2] = -V3[1].imag();
  P3[3] = -V3[0].imag();
  denom = COUP / ((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) -
                  (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-cI) *
          (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]);
  V3[3] = denom * (-cI) *
          (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3] * F2[4]);
  V3[4] = denom * (-cI) *
          (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) +
           cI * (F1[3] * F2[4] + F1[4] * F2[3]));
  V3[5] = denom * (-cI) *
          (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5] * F2[3]);
}

void FFV4_0(extras::complex F1[], extras::complex F2[],
            extras::complex V3[], extras::complex COUP,
            extras::complex &vertex) {
  static extras::complex cI = extras::complex(0., 1.);
  extras::complex TMP0;
  extras::complex TMP3;
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
          F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  TMP3 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
          F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  vertex = COUP * (-1.) * (+cI * (TMP0) + 2. * cI * (TMP3));
}

void FFV4_3(extras::complex F1[], extras::complex F2[],
            extras::complex COUP, double M3, double W3,
            extras::complex V3[]) {
  static extras::complex cI = extras::complex(0., 1.);
  extras::complex denom;
  extras::complex TMP1;
  double P3[4];
  double OM3;
  extras::complex TMP4;
  OM3 = 0.;
  if (M3 != 0.)
    OM3 = 1. / (M3 * M3);
  V3[0] = +F1[0] + F2[0];
  V3[1] = +F1[1] + F2[1];
  P3[0] = -V3[0].real();
  P3[1] = -V3[1].real();
  P3[2] = -V3[1].imag();
  P3[3] = -V3[0].imag();
  TMP4 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
          F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
          F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  denom = COUP / ((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) -
                  (P3[3] * P3[3]) - M3 * (M3 - cI * W3));
  V3[2] = denom * (-2. * cI) *
          (OM3 * -1. / 2. * P3[0] * (TMP1 + 2. * (TMP4)) +
           (+1. / 2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] +
            F1[5] * F2[3]));
  V3[3] = denom * (-2. * cI) *
          (OM3 * -1. / 2. * P3[1] * (TMP1 + 2. * (TMP4)) +
           (-1. / 2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] +
            F1[5] * F2[2]));
  V3[4] = denom * 2. * cI *
          (OM3 * 1. / 2. * P3[2] * (TMP1 + 2. * (TMP4)) +
           (+1. / 2. * cI * (F1[2] * F2[5]) - 1. / 2. * cI * (F1[3] * F2[4]) -
            cI * (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
  V3[5] = denom * 2. * cI *
          (OM3 * 1. / 2. * P3[3] * (TMP1 + 2. * (TMP4)) +
           (+1. / 2. * (F1[2] * F2[4]) - 1. / 2. * (F1[3] * F2[5]) -
            F1[4] * F2[2] + F1[5] * F2[3]));
}

void FFV2_0(extras::complex F1[], extras::complex F2[],
            extras::complex V3[], extras::complex COUP,
            extras::complex &vertex) {
  static extras::complex cI = extras::complex(0., 1.);
  extras::complex TMP0;
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
          F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * -cI * TMP0;
}

void FFV2_4_0(extras::complex F1[], extras::complex F2[],
              extras::complex V3[], extras::complex COUP1,
              extras::complex COUP2, extras::complex &vertex) {
  extras::complex tmp;
  FFV2_0(F1, F2, V3, COUP1, vertex);
  FFV4_0(F1, F2, V3, COUP2, tmp);
  vertex = vertex + tmp;
}

} // namespace MG5_sm
