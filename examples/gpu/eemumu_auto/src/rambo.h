#include <vector>
#include "mg5Complex.h"

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

std::vector<std::vector<floa_t *>> get_momenta(int ninitial, double energy,
                                               std::vector<floa_t> masses,
                                               double &wgt, int dim);

std::vector<floa_t *> rambo(double et, std::vector<floa_t> &xm, double &wt);
