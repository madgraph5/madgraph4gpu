#include <CL/sycl.hpp>
#include <cmath>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>

#include "rambo.h"

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
  if (iranmr < 0)
    iranmr = 97;
  if (jranmr < 0)
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


bool pass_cuts(const std::vector<double*>& pout){
  // return true if the events pass the cuts. (fails if not)
  // two cut are implemented
  // - pt cut on gluon (particle 3 and 4): 20 GeV
  // - delta_R between the gluon (0.1)

  // pt cut on first gluon
  double pt2 =pout[2][1]*pout[2][1] + pout[2][2]*pout[2][2];
  if (pt2 < 400){
    return false;
  }
  
  // pt cut on second gluon
  pt2 =pout[3][1]*pout[3][1] + pout[3][2]*pout[3][2];
  if (pt2 < 400){
    return false;
  }  
   
  //delta_R between the gluon
  // delta_R**2 = delta_phi**2 + delta_eta**2
  // computing delta_phi first
  double denom = sqrt(pout[2][1]*pout[2][1] + pout[2][2]*pout[2][2]) * sqrt(pout[3][1]*pout[3][1] + pout[3][2]*pout[3][2]);
  double delta_phi = acos((pout[2][1]*pout[3][1] + pout[2][2]*pout[3][2]) / denom);
  
  // computing delta_eta (note rap1/2 are not lorentz invariant but the difference is.
  // This evaluates rapidity in the center of mass frame (not the same as lab frame)
  double rap1 = 0.5 * log( (pout[2][0] + pout[2][3])/(pout[2][0] - pout[2][3]));
  double rap2 = 0.5 * log( (pout[3][0] + pout[3][3])/(pout[3][0] - pout[3][3]));

  //computing deltaR and check bounds
  double delta_rap2 = delta_phi*delta_phi + (rap1 -rap2) * (rap1 - rap2);
  if (delta_rap2 < 0.01){
    return false;
  }
    
  
  return true;
}

std::vector<std::vector<double *>> get_momenta(int ninitial, double energy,
                                               std::vector<double> masses,
                                               double &wgt, int dim) {
  //---- auxiliary function to change convention between MadGraph5_aMC@NLO and
  // rambo
  //---- four momenta.
  int nexternal = masses.size();
  int nfinal = nexternal - ninitial;
  double e2 = pow(energy, 2);
  double m1 = masses[0];
  std::vector<std::vector<double *>> p2;

  if (ninitial == 1) {
	  int done = 0;
      while (done < dim){  
			 //    for (int d = 0; d < dim; ++d) {
      // Momenta for the incoming particle
      std::vector<double *> p(1, new double[4]);
      p[0][0] = m1;
      p[0][1] = 0.;
      p[0][2] = 0.;
      p[0][3] = 0.;

      std::vector<double> finalmasses(++masses.begin(), masses.end());
      std::vector<double *> p_rambo = rambo(m1, finalmasses, wgt);
      if (pass_cuts(p_rambo)){
	p.insert(++p.begin(), p_rambo.begin(), p_rambo.end());
	p2.push_back(p);
	done++; 
      }
    }
    return p2;
  }

  else if (ninitial != 2) {
    std::cout << "Rambo needs 1 or 2 incoming particles" << std::endl;
    exit(-1);
  }

  if (nfinal == 1)
    energy = m1;

  double m2 = masses[1];

  double mom =
      sqrt((pow(e2, 2) - 2 * e2 * pow(m1, 2) + pow(m1, 4) -
            2 * e2 * pow(m2, 2) - 2 * pow(m1, 2) * pow(m2, 2) + pow(m2, 4)) /
           (4 * e2));
  double energy1 = sqrt(pow(mom, 2) + pow(m1, 2));
  double energy2 = sqrt(pow(mom, 2) + pow(m2, 2));
  // Set momenta for incoming particles
  for (int d = 0; d < dim; ++d) {

    std::vector<double *> p(1, new double[4]);
    p[0][0] = energy1;
    p[0][1] = 0;
    p[0][2] = 0;
    p[0][3] = mom;
    p.push_back(new double[4]);
    p[1][0] = energy2;
    p[1][1] = 0;
    p[1][2] = 0;
    p[1][3] = -mom;

    if (nfinal == 1) {
      p.push_back(new double[4]);
      p[2][0] = energy;
      wgt = 1;
      p2.push_back(p);
    }
    std::vector<double> finalmasses(++(++masses.begin()), masses.end());
    std::vector<double *> p_rambo = rambo(energy, finalmasses, wgt);
    p.insert(++(++p.begin()), p_rambo.begin(), p_rambo.end());
    p2.push_back(p);
  }
  return p2;
}

std::vector<double *> rambo(double et, std::vector<double> &xm, double &wt) {
  /**********************************************************************
   *                       rambo                                         *
   *    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)                  *
   *                                                                     *
   *    a democratic multi-particle phase space generator                *
   *    authors:  s.d. ellis,  r. kleiss,  w.j. stirling                 *
   *    this is version 1.0 -  written by r. kleiss                      *
   *    -- adjusted by hans kuijf, weights are logarithmic (20-08-90)    *
   *                                                                     *
   *    n  = number of particles                                         *
   *    et = total centre-of-mass energy                                 *
   *    xm = particle masses ( dim=nexternal-nincoming )                 *
   *    p  = particle momenta ( dim=(4,nexternal-nincoming) )            *
   *    wt = weight of the event                                         *
   ***********************************************************************/
  int n = xm.size();
  std::vector<double *> q, p;
  std::vector<double> z(n), r(4), b(3), p2(n), xm2(n), e(n), v(n);
  static std::vector<int> iwarn(5, 0);
  static double acc = 1e-14;
  static int itmax = 6, ibegin = 0;
  static double twopi = 8. * atan(1.);
  static double po2log = log(twopi / 4.);

  for (int i = 0; i < n; i++) {
    q.push_back(new double[4]);
    p.push_back(new double[4]);
  }
  // initialization step: factorials for the phase space weight
  if (ibegin == 0) {
    ibegin = 1;
    z[1] = po2log;
    for (int k = 2; k < n; k++)
      z[k] = z[k - 1] + po2log - 2. * log(double(k - 1));
    for (int k = 2; k < n; k++)
      z[k] = (z[k] - log(double(k)));
  }
  // check on the number of particles
  if (n < 1 || n > 101) {
    std::cout << "Too few or many particles: " << n << std::endl;
    exit(-1);
  }
  // check whether total energy is sufficient; count nonzero masses
  double xmt = 0.;
  int nm = 0;
  for (int i = 0; i < n; i++) {
    if (xm[i] != 0.)
      nm = nm + 1;
    xmt = xmt + sycl::abs(xm[i]);
  }
  if (xmt > et) {
    std::cout << "Too low energy: " << et << " needed " << xmt << std::endl;
    exit(-1);
  }
  // the parameter values are now accepted

  // generate n massless momenta in infinite phase space
  for (int i = 0; i < n; i++) {
    double r1 = rn(1);
    double c = 2. * r1 - 1.;
    double s = sqrt(1. - c * c);
    double f = twopi * rn(2);
    r1 = rn(3);
    double r2 = rn(4);
    q[i][0] = -log(r1 * r2);
    q[i][3] = q[i][0] * c;
    q[i][2] = q[i][0] * s * cos(f);
    q[i][1] = q[i][0] * s * sin(f);
  }
  // calculate the parameters of the conformal transformation
  for (int i = 0; i < 4; i++)
    r[i] = 0.;
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 4; k++)
      r[k] = r[k] + q[i][k];
  }
  double rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
  for (int k = 1; k < 4; k++)
    b[k - 1] = -r[k] / rmas;
  double g = r[0] / rmas;
  double a = 1. / (1. + g);
  double x = et / rmas;

  // transform the q's conformally into the p's
  for (int i = 0; i < n; i++) {
    double bq = b[0] * q[i][1] + b[1] * q[i][2] + b[2] * q[i][3];
    for (int k = 1; k < 4; k++)
      p[i][k] = x * (q[i][k] + b[k - 1] * (q[i][0] + a * bq));
    p[i][0] = x * (g * q[i][0] + bq);
  }

  for (int i = 0; i < n; ++i) {
    delete[] q[i];
  }

  // calculate weight and possible warnings
  wt = po2log;
  if (n != 2)
    wt = (2. * n - 4.) * log(et) + z[n - 1];
  if (wt < -180.) {
    if (iwarn[0] <= 5)
      std::cout << "Too small wt, risk for underflow: " << wt << std::endl;
    iwarn[0] = iwarn[0] + 1;
  }
  if (wt > 174.) {
    if (iwarn[1] <= 5)
      std::cout << "Too large wt, risk for overflow: " << wt << std::endl;
    iwarn[1] = iwarn[1] + 1;
  }

  // return for weighted massless momenta
  if (nm == 0) {
    // return log of weight
    return p;
  }

  // massive particles: rescale the momenta by a factor x
  double xmax = sqrt(1. - pow(xmt / et, 2));
  for (int i = 0; i < n; i++) {
    xm2[i] = pow(xm[i], 2);
    p2[i] = pow(p[i][0], 2);
  }
  int iter = 0;
  x = xmax;
  double accu = et * acc;
  while (true) {
    double f0 = -et;
    double g0 = 0.;
    double x2 = x * x;
    for (int i = 0; i < n; i++) {
      e[i] = sqrt(xm2[i] + x2 * p2[i]);
      f0 = f0 + e[i];
      g0 = g0 + p2[i] / e[i];
    }
    if (sycl::abs(f0) <= accu)
      break;
    iter = iter + 1;
    if (iter > itmax) {
      std::cout << "Too many iterations without desired accuracy: " << itmax
                << std::endl;
      break;
    }
    x = x - f0 / (x * g0);
  }
  for (int i = 0; i < n; i++) {
    v[i] = x * p[i][0];
    for (int k = 1; k < 4; k++)
      p[i][k] = x * p[i][k];
    p[i][0] = e[i];
  }

  // calculate the mass-effect weight factor
  double wt2 = 1.;
  double wt3 = 0.;
  for (int i = 0; i < n; i++) {
    wt2 = wt2 * v[i] / e[i];
    wt3 = wt3 + pow(v[i], 2) / e[i];
  }
  double wtm = (2. * n - 3.) * log(x) + log(wt2 / wt3 * et);

  // return for  weighted massive momenta
  wt = wt + wtm;
  if (wt < -180.) {
    if (iwarn[2] <= 5)
      std::cout << "Too small wt, risk for underflow: " << wt << std::endl;
    iwarn[2] = iwarn[2] + 1;
  }
  if (wt > 174.) {
    if (iwarn[3] <= 5)
      std::cout << "Too large wt, risk for overflow: " << wt << std::endl;
    iwarn[3] = iwarn[3] + 1;
  }
  // return log of weight
  return p;
}
