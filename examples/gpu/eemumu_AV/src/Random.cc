#include "Random.h"

double Random::ranmar() {
  /*     -----------------
   * universal random number generator proposed by marsaglia and zaman
   * in report fsu-scri-87-50
   * in this version rvec is a double precision variable. */
  double uni = ranu[iranmr] - ranu[jranmr];
  if (uni < 0)
    uni = uni + 1;
  ranu[iranmr] = uni;
  iranmr = iranmr - 1;
  jranmr = jranmr - 1;
  if (iranmr == 0)
    iranmr = 97;
  if (jranmr == 0)
    jranmr = 97;
  ranc = ranc - rancd;
  if (ranc < 0)
    ranc = ranc + rancm;
  uni = uni - ranc;
  if (uni < 0)
    uni = uni + 1;
  return uni;
}

void Random::rmarin(int ij, int kl) {
  /*     -----------------
   * initializing routine for ranmar, must be called before generating
   * any pseudorandom numbers with ranmar. the input values should be in
   * the ranges 0<=ij<=31328 ; 0<=kl<=30081 */
  /* this shows correspondence between the simplified input seeds ij, kl
   * and the original marsaglia-zaman seeds i,j,k,l.
   * to get the standard values in the marsaglia-zaman paper (i=12,j=34
   * k=56,l=78) put ij=1802, kl=9373 */
  int i = ij / 177 % 177 + 2;
  int j = ij % 177 + 2;
  int k = (kl / 169) % 178 + 1;
  int l = kl % 169;
  for (int ii = 1; ii < 98; ii++) {
    double s = 0;
    double t = .5;
    for (int jj = 1; jj < 25; jj++) {
      int m = ((i * j % 179) * k) % 179;
      i = j;
      j = k;
      k = m;
      l = (53 * l + 1) % 169;
      if ((l * m) % 64 >= 32)
        s = s + t;
      t = .5 * t;
    }
    ranu[ii] = s;
  }
  ranc = 362436. / 16777216.;
  rancd = 7654321. / 16777216.;
  rancm = 16777213. / 16777216.;
  iranmr = 97;
  jranmr = 33;
}

double rn(int idummy) {
  static Random rand;
  double ran;
  static int init = 1;
  // Prevent unused variable warning
  if (false)
    idummy = idummy;
  if (init == 1) {
    init = 0;
    rand.rmarin(1802, 9373);
  }

  while (true) {
    ran = rand.ranmar();
    if (ran > 1e-16)
      break;
  }
  return ran;
}
