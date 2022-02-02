#include "bridge.h"

#ifdef __CUDACC__

/**
 * The C++ function instantiating the Bridge class with all necessary parameters
 * and calling the host_sequence function to execute the Cuda/C++ code for the
 * matrix element calculation
 */
void cubridge(double *momenta, double *mes) {
  Bridge<double> b = Bridge<double>(16, 4, 4, 4, 16);
  b.gpu_sequence(momenta, mes);
}

extern "C" {
/**
 * The C symbol being called from the fortran code (in audo_dsig1.f)
 */
void fcubridge_(double *mom, double *mes) { cubridge(mom, mes); }
}

#else

/**
 * The C++ function instantiating the Bridge class with all necessary parameters
 * and calling the host_sequence function to execute the Cuda/C++ code for the
 * matrix element calculation
 */
void cppbridge(double *momenta, double *mes) {
  Bridge<double> b = Bridge<double>(16, 4, 4, 4, 16);
  b.cpu_sequence(momenta, mes);
}

extern "C" {
/**
 * The C symbol being called from the fortran code (in audo_dsig1.f)
 */
void fcppbridge_(double *mom, double *mes) { cppbridge(mom, mes); }
}

#endif // __CUDACC__
