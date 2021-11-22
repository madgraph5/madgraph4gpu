#include "transpose.h"
#include <iostream>

// m == DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE)
// NEXTERNAL = 4 (nexternal.inc)
// NB_PAGE = 16 (vector.inc)
// ncomb = 16 (mgOnGpuConfig.h)

struct CudaInit {
  CudaInit();
};

CudaInit::CudaInit() { std::cout << "cuda init" << std::endl; }

static CudaInit cuInit;

void bridge(double *momenta, double *mes) {

  // par 1: number of events (NB_PAGE)
  // par 2: number of particles / event (NEXTERNAL)
  // par 3: number of momenta / particle
  // par 4: stride length
  // par 5: number of good helicities (ncomb)
  Matrix<double> t = Matrix<double>(16, 4, 4, 4, 16);
  t.hst_transpose(momenta, mes);
}

extern "C" {
void bridge_(double *momenta, double *mes) {
  bridge(momenta, mes);
  // for (int i = 0; i < 16; ++i) {
  //   printf("%f x ", mes[i]);
  // }
}
}
