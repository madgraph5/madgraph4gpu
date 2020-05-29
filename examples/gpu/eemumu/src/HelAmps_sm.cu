#ifndef gHelAmps_sm
#define gHelAmps_sm

#include "HelAmps_sm.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include <cuda_runtime.h>
#include <thrust/complex.h>

#define gpuErrchk2(ans)                                                        \
  { gpuAssert2((ans), __FILE__, __LINE__); }

__device__ void gpuAssert2(cudaError_t code, const char *file, int line,
                           bool abort = true) {
  if (code != cudaSuccess) {
    printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
  }
}

extern __constant__ int cHel[16][4];
extern __constant__ double cMME[4];
extern __constant__ int cPerm[4];
extern __constant__ thrust::complex<double> cIPC[3];
extern __constant__ double cIPD[2];
// extern __shared__ thrust::complex<double> sw[6][6];

namespace gMG5_sm {

#ifdef DEBUG
__device__ void debugMsg(const char *msg) {
  printf("%i.%i-%s ", blockIdx.x, threadIdx.x, msg);
}
#endif

// later
/*
void CPPProcess::sigmaKin(int ncomb, bool (&goodhel)[16], int &ntry,
                          int &sum_hel, int &ngood, int (&igood)[16],
                          int &jhel) {
                    */
__global__ void sigmaKin(cudaPitchedPtr tp, double *meDevPtr, size_t mePitch,
                         bool debug, bool verbose) {

  // for (int xx = 0; xx < 384; ++xx) {

  int nprocesses = 1;

  int dim = blockIdx.x * blockDim.x + threadIdx.x;

  char *devPtr = (char *)tp.ptr;
  size_t dpt = tp.pitch;
  size_t slicePitch = dpt * 4;

  char *dps = devPtr + dim * slicePitch;
  double *matrix_element = (double *)((char *)meDevPtr + dim * mePitch);

  thrust::complex<double> amp[2];
  double t[1];

  // <later>
  // Local variables and constants
  const int ncomb = 16;
  static bool goodhel[ncomb] = {ncomb * false};
  static int ntry = 0, ngood = 0;
  static int igood[ncomb];
  // </later>

  // Denominators: spins, colors and identical particles
  const int denominators[1] = {4}; // nprocesses

  ntry = ntry + 1;

  // Reset the matrix elements
  for (int i = 0; i < nprocesses; ++i) { // nprocesses
    matrix_element[i] = 0.;
  }

  // sr fixme // better to run the first n calculations serial?
  // if (sum_hel == 0 || ntry < 10) {
  // Calculate the matrix element for all helicities

  for (int ihel = 0; ihel < ncomb; ihel++) {
    if (goodhel[ihel] || ntry < 2) {

      calculate_wavefunctions(ihel, dps, dpt, amp, debug, verbose);
      matrix_1_epem_mupmum(t[0], amp);

      double tsum = 0;
      for (int iproc = 0; iproc < nprocesses; iproc++) {
        matrix_element[iproc] += t[iproc];
        tsum += t[iproc];
      }

      // Store which helicities give non-zero result
      if (tsum != 0. && !goodhel[ihel]) {
        goodhel[ihel] = true;
        ngood++;
        igood[ngood] = ihel;
      }
    }
  }

  for (int i = 0; i < nprocesses; ++i) {
    matrix_element[i] /= denominators[i];
  }
  // }

  // printf("%d - %e\n", dim, t[0]);
}

// --> calculate multi-dimensional amp
__device__ void matrix_1_epem_mupmum(double &matrix,
                                     thrust::complex<double> amp[2]) {
  int i, j;
  // Local variables
  // const int ngraphs = 2;
  const int ncolor = 1;
  thrust::complex<double> ztemp;
  thrust::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{1}};

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1];

  // Sum and square the color flows to get the matrix element
  matrix = 0;
  for (i = 0; i < ncolor; i++) {
    ztemp = 0.;
    for (j = 0; j < ncolor; j++)
      ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + (ztemp * conj(jamp[i])).real() / denom[i];
  }

  // Store the leading color flows for choice of color
  // sr fixme // maybe this needs to go outside the loop? does it need a
  // dimension?
  /* sr fixme
  for (i = 0; i < ncolor; i++)
    jamp2[0][i] += (jamp[i] * conj(jamp[i])).real();
  */
}

__device__ void calculate_wavefunctions(int ihel, char *dps, size_t dpt,
                                        thrust::complex<double> amp[2],
                                        bool debug, bool verbose) {
#ifdef DEBUG
  debugMsg("%>");
  if (debug) {
    printf("\n\nblock (%i / %i), thread (%i)\n\n", blockIdx.x, blockDim.x,
           threadIdx.x);
  }
#endif

  double ZERO = 0.00;
  thrust::complex<double> sw[6][6];

  // Calculate all wavefunctions
  oxxxxx((double *)(dps + cPerm[0] * dpt), cMME[0], cHel[ihel][0], -1, sw[0]);
  ixxxxx((double *)(dps + cPerm[1] * dpt), cMME[1], cHel[ihel][1], +1, sw[1]);
  ixxxxx((double *)(dps + cPerm[2] * dpt), cMME[2], cHel[ihel][2], -1, sw[2]);
  oxxxxx((double *)(dps + cPerm[3] * dpt), cMME[3], cHel[ihel][3], +1, sw[3]);
  FFV1P0_3(sw[1], sw[0], cIPC[0], ZERO, ZERO, sw[4]);
  FFV2_4_3(sw[1], sw[0], -cIPC[1], cIPC[2], cIPD[0], cIPD[1], sw[5]);
  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(sw[2], sw[3], sw[4], cIPC[0], &amp[0]);
  FFV2_4_0(sw[2], sw[3], sw[5], -cIPC[1], cIPC[2], &amp[1]);

#ifdef DEBUG
  if (debug) {
    printf("\n\n >>> DEBUG >>> DEBUG >>> DEBUG >>>\n");

    printf("\nHelicities: %d %d %d %d\n", cHel[ihel][0], cHel[ihel][1],
           cHel[ihel][2], cHel[ihel][3]);

    printf("\nMomenta:\n");
    for (int i = 0; i < 4; ++i) {
      printf("%i %e %e %e %e\n", i, dp[cPerm[i]][0], dp[cPerm[i]][1],
             dp[cPerm[i]][2], dp[cPerm[i]][3]);
    }

    printf("\nMasses: %e, %e, %e, %e\n", cMME[0], cMME[1], cMME[2], cMME[3]);

    printf("\nAmplitudes: (%e, %e), (%e, %e)\n", amp[0].real(), amp[0].imag(),
           amp[1].real(), amp[1].imag());

    printf("\nWavefuncs:\n");
    for (int i = 0; i < 6; ++i) {
      printf("%i ", i);
      for (int j = 0; j < 6; ++j) {
        double re = sw[i][j].real(), im = sw[i][j].imag();
        if (re == 0 && im == 0) {
          printf("0, ");
        } else {
          printf("(%e, %e), ", re, im);
        }
      }
      printf("\n");
    }

    printf("\n\n <<< DEBUG <<< DEBUG <<< DEBUG <<<\n\n");
  }
  debugMsg("<%");
#endif
}

__device__ void ixxxxx(double p[4], double fmass, int nhel, int nsf,
                       thrust::complex<double> fi[6]) {
#ifdef DEBUG
  debugMsg("b>");
#endif
  thrust::complex<double> chi[2];
  double sqp0p3;
  int nh;
  fi[0] = thrust::complex<double>(-p[0] * nsf, -p[3] * nsf);
  fi[1] = thrust::complex<double>(-p[1] * nsf, -p[2] * nsf);
  nh = nhel * nsf;
  if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0) {
    sqp0p3 = 0.0;
  } else {
    sqp0p3 = sqrt(max(p[0] + p[3], 0.0)) * nsf;
  }
  chi[0] = thrust::complex<double>(sqp0p3, 0.0);
  if (sqp0p3 == 0.0) {
    chi[1] = thrust::complex<double>(-nhel * sqrt(2.0 * p[0]), 0.0);
  } else {
    chi[1] = thrust::complex<double>(nh * p[1], p[2]) / sqp0p3;
  }
  if (nh == 1) {
    fi[2] = thrust::complex<double>(0.0, 0.0);
    fi[3] = thrust::complex<double>(0.0, 0.0);
    fi[4] = chi[0];
    fi[5] = chi[1];
  } else {
    fi[2] = chi[1];
    fi[3] = chi[0];
    fi[4] = thrust::complex<double>(0.0, 0.0);
    fi[5] = thrust::complex<double>(0.0, 0.0);
  }
#ifdef DEBUG
  debugMsg("<b");
#endif
  return;
}

__device__ void oxxxxx(double p[4], double fmass, int nhel, int nsf,
                       thrust::complex<double> fo[6]) {
#ifdef DEBUG
  debugMsg("a>");
#endif
  thrust::complex<double> chi[2];
  double sqp0p3;
  int nh;
  fo[0] = thrust::complex<double>(p[0] * nsf, p[3] * nsf);
  fo[1] = thrust::complex<double>(p[1] * nsf, p[2] * nsf);
  nh = nhel * nsf;
  if ((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00)) {
    sqp0p3 = 0.00;
  } else {
    sqp0p3 = sqrt(max(p[0] + p[3], 0.00)) * nsf;
  }
  chi[0] = thrust::complex<double>(sqp0p3, 0.00);
  if (sqp0p3 == 0.000) {
    chi[1] = thrust::complex<double>(-nhel, 0.00) * sqrt(2.0 * p[0]);
  } else {
    chi[1] = thrust::complex<double>(nh * p[1], -p[2]) / sqp0p3;
  }
  if (nh == 1) {
    fo[2] = chi[0];
    fo[3] = chi[1];
    fo[4] = thrust::complex<double>(0.00, 0.00);
    fo[5] = thrust::complex<double>(0.00, 0.00);
  } else {
    fo[2] = thrust::complex<double>(0.00, 0.00);
    fo[3] = thrust::complex<double>(0.00, 0.00);
    fo[4] = chi[1];
    fo[5] = chi[0];
  }
#ifdef DEBUG
  debugMsg("<a");
#endif
  return;
}

__device__ void FFV1_0(thrust::complex<double> F1[],
                       thrust::complex<double> F2[],
                       thrust::complex<double> V3[],
                       thrust::complex<double> COUP,
                       thrust::complex<double> *vertex) {
#ifdef DEBUG
  debugMsg("g>");
#endif
  thrust::complex<double> cI = thrust::complex<double>(0., 1.);
  thrust::complex<double> TMP2 =
      (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
       (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
        (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
         F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  (*vertex) = COUP * -cI * TMP2;
#ifdef DEBUG
  debugMsg("<g");
#endif
}

__device__ void FFV2_3(thrust::complex<double> F1[],
                       thrust::complex<double> F2[],
                       thrust::complex<double> COUP, double M3, double W3,
                       thrust::complex<double> V3[]) {
#ifdef DEBUG
  debugMsg("e>");
#endif
  thrust::complex<double> cI = thrust::complex<double>(0., 1.);
  thrust::complex<double> denom;
  thrust::complex<double> TMP1;
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
#ifdef DEBUG
  debugMsg("<e");
#endif
}

__device__ void FFV4_3(thrust::complex<double> F1[],
                       thrust::complex<double> F2[],
                       thrust::complex<double> COUP, double M3, double W3,
                       thrust::complex<double> V3[]) {
#ifdef DEBUG
  debugMsg("f>");
#endif
  thrust::complex<double> cI = thrust::complex<double>(0., 1.);
  thrust::complex<double> denom;
  thrust::complex<double> TMP1;
  double P3[4];
  double OM3;
  thrust::complex<double> TMP4;
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
#ifdef DEBUG
  debugMsg("<f");
#endif
}

__device__ void FFV2_4_3(thrust::complex<double> F1[],
                         thrust::complex<double> F2[],
                         thrust::complex<double> COUP1,
                         thrust::complex<double> COUP2, double M3, double W3,
                         thrust::complex<double> V3[]) {
#ifdef DEBUG
  debugMsg("d>");
#endif
  int i;
  thrust::complex<double> *Vtmp;
  gpuErrchk2(cudaMalloc(&Vtmp, 6 * sizeof(thrust::complex<double>)));
  *Vtmp = thrust::complex<double>(0, 0);
  FFV2_3(F1, F2, COUP1, M3, W3, V3);
  FFV4_3(F1, F2, COUP2, M3, W3, Vtmp);
  i = 2;
  while (i < 6) {
    V3[i] = V3[i] + Vtmp[i];
    i++;
  }
  gpuErrchk2(cudaFree(Vtmp));
#ifdef DEBUG
  debugMsg("<d");
#endif
}

__device__ void FFV1P0_3(thrust::complex<double> F1[],
                         thrust::complex<double> F2[],
                         thrust::complex<double> COUP, double M3, double W3,
                         thrust::complex<double> V3[]) {
#ifdef DEBUG
  debugMsg("c>");
#endif
  thrust::complex<double> cI = thrust::complex<double>(0., 1.);
  double P3[4];
  thrust::complex<double> denom;
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
#ifdef DEBUG
  debugMsg("<c");
#endif
}

__device__ void FFV4_0(thrust::complex<double> F1[],
                       thrust::complex<double> F2[],
                       thrust::complex<double> V3[],
                       thrust::complex<double> COUP,
                       thrust::complex<double> *vertex) {
#ifdef DEBUG
  debugMsg("j>");
#endif
  thrust::complex<double> cI = thrust::complex<double>(0., 1.);
  thrust::complex<double> TMP0, TMP3;
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
          F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  TMP3 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
          F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  (*vertex) = COUP * (-1.) * (+cI * (TMP0) + 2. * cI * (TMP3));
#ifdef DEBUG
  debugMsg("<j");
#endif
}

__device__ void FFV2_0(thrust::complex<double> F1[],
                       thrust::complex<double> F2[],
                       thrust::complex<double> V3[],
                       thrust::complex<double> COUP,
                       thrust::complex<double> *vertex) {
#ifdef DEBUG
  debugMsg("i>");
#endif
  thrust::complex<double> cI = thrust::complex<double>(0., 1.);
  thrust::complex<double> TMP0;
  TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
          F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  (*vertex) = COUP * -cI * TMP0;
#ifdef DEBUG
  debugMsg("<i");
#endif
}

__device__ void
FFV2_4_0(thrust::complex<double> F1[], thrust::complex<double> F2[],
         thrust::complex<double> V3[], thrust::complex<double> COUP1,
         thrust::complex<double> COUP2, thrust::complex<double> *vertex) {
#ifdef DEBUG
  debugMsg("h>");
#endif
  thrust::complex<double> tmp;
  FFV2_0(F1, F2, V3, COUP1, vertex);
  FFV4_0(F1, F2, V3, COUP2, &tmp);
  (*vertex) = (*vertex) + tmp;
#ifdef DEBUG
  debugMsg("<h");
#endif
}

} // namespace gMG5_sm

#endif // gHelAmps_sm
