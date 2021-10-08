#include <iostream>
#include "transpose.h"

// m == DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE)
// NEXTERNAL = 4 (nexternal.inc)
// NB_PAGE = 16 (vector.inc)

void bridge(double *m) {

  std::cout << "kernel,";
  Matrix<double> t(16, 4, 4, 2);
  t.hst_transpose(m);

}

extern "C" {

  void bridge_(double *m) {
    bridge(m);
  }

}
