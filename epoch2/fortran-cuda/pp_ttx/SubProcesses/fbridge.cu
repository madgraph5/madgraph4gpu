#include "bridge.h"

/**
 * The C++ function instantiating the Bridge class with all necessary parameters
 * and calling the host_sequence function to execute the Cuda/C++ code for the
 * matrix element calculation
 */
void fbridge(double *momenta, double *mes) {
  Bridge<double> b = Bridge<double>(16, 4, 4, 4, 16);
  b.host_sequence(momenta, mes);
}

extern "C" {
/**
 * The C symbol being called from the fortran code (in audo_dsig1.f)
 */
void fbridge_(double *momenta, double *mes) { fbridge(momenta, mes); }
}
