#include <iostream>

// m == DOUBLE PRECISION P_MULTI(0:3, NEXTERNAL, NB_PAGE)
// NEXTERNAL = 4 (nexternal.inc)
// NB_PAGE = 16 (vector.inc)

void bridge(double *m) {

  std::cout << m[0] << ",";

}

extern "C" {

  void bridge_(double *m) {
    bridge(m);
  }

}
