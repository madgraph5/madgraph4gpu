#pragma once
#include "extras.h"

class Random {
public:
  double ranmar();
  void rmarin(int ij, int kl);

private:
  double ranu[98];
  double ranc, rancd, rancm;
  int iranmr, jranmr;
};

double rn(int idummy);

double** get_momenta(int ninitial, double energy, double masses[],
                             double &wgt);

double** rambo(double et, double xm[2], double &wt);
