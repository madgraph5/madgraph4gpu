#include "bridge.h"

void fbridge(double *momenta, double *mes) {
  Bridge<double> b = Bridge<double>(16, 4, 4, 4, 16);
  b.host_sequence(momenta, mes);
}

extern "C" {
void fbridge_(double *momenta, double *mes) { fbridge(momenta, mes); }
}
