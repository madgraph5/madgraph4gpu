#include <cmath>
#include <cstdlib>
#include <iostream>

#include "vrambo.h"
#include "Random.h"

//std::vector<std::vector<double *>> // output is an AOS: momenta[nevt][nexternal][4]
//get_momenta(int ninitial, double energy, const std::vector<double> masses, double &wgt, int nevt) {

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta.
void get_momenta( const int ninitial,    // input: #particles_initial
                  const double energy,   // input: energy
                  const double masses[], // input: masses[npar]
                  double momenta1d[],    // output: momenta[nevt][npar][4] as an AOS
                  double wgts[],         // output: wgts[nevt]
                  const int npar,        // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt )       // input: #events
{
  const int nexternal = npar;
  const int nfinal = nexternal - ninitial;
  const double e2 = pow(energy, 2);
  const double m1 = masses[0];

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)

  double (*momenta)[npar][np4] = (double (*)[npar][np4]) momenta1d; // cast to multiD array pointer

  // #Initial==1
  if (ninitial == 1) {
    for (int ievt = 0; ievt < nevt; ++ievt) {
      // Momenta for the incoming particle
      momenta[ievt][0][0] = m1;
      momenta[ievt][0][1] = 0;
      momenta[ievt][0][2] = 0;
      momenta[ievt][0][3] = 0;
      // Momenta for the outgoing particles
      std::vector<double> finalmasses;
      for (int ipar = 0; ipar < npar-ninitial; ++ipar) // loop over npar-ninitial particles from rambo
        finalmasses.push_back( masses[ipar] );
      double wgt;
      std::vector<double *> p_rambo = rambo(m1, finalmasses, wgt);
      for (int ipar = 0; ipar < npar-ninitial; ++ipar) // loop over npar-ninitial particles from rambo
        for (int ip4 = 0; ip4 < np4; ++ip4)
          momenta[ievt][ipar+ninitial][ip4] = p_rambo[ipar][ip4];
      // Event weight
      wgts[ievt] = wgt;
      // Match the 'new double[4]' in rambo
      for (int ipar = 0; ipar < npar-ninitial; ++ipar) // loop over npar-ninitial particles from rambo
        delete[] p_rambo[ipar];
    }
    return;
  }

  // #Initial>2 (error)
  else if (ninitial != 2) {
    std::cout << "Rambo needs 1 or 2 incoming particles" << std::endl;
    exit(-1);
  }

  // #Initial==2
  else {
    double m2 = masses[1];
    double mom =
      sqrt((pow(e2, 2) - 2 * e2 * pow(m1, 2) + pow(m1, 4) -
            2 * e2 * pow(m2, 2) - 2 * pow(m1, 2) * pow(m2, 2) + pow(m2, 4)) /
           (4 * e2));
    double energy1 = sqrt(pow(mom, 2) + pow(m1, 2));
    double energy2 = sqrt(pow(mom, 2) + pow(m2, 2));
    for (int ievt = 0; ievt < nevt; ++ievt) {
      // Momenta for the incoming particles
      momenta[ievt][0][0] = energy1;
      momenta[ievt][0][1] = 0;
      momenta[ievt][0][2] = 0;
      momenta[ievt][0][3] = mom;
      momenta[ievt][1][0] = energy2;
      momenta[ievt][1][1] = 0;
      momenta[ievt][1][2] = 0;
      momenta[ievt][1][3] = -mom;
      // #Initial==2, #Final==1
      if (nfinal == 1) {
        // Momenta for the outgoing particle
        momenta[ievt][2][0] = m1;
        momenta[ievt][2][1] = 0;
        momenta[ievt][2][2] = 0;
        momenta[ievt][2][3] = 0;
        // Event weight
        double wgt = 1;
        wgts[ievt] = wgt;
      }
      // #Initial==2, #Final>1
      else {
        // Momenta for the outgoing particles
        std::vector<double> finalmasses;
        for (int ipar = 0; ipar < npar-ninitial; ++ipar) // loop over npar-ninitial particles from rambo
          finalmasses.push_back( masses[ipar] );
        double wgt;
        std::vector<double *> p_rambo = rambo(energy, finalmasses, wgt);
        for (int ipar = 0; ipar < npar-ninitial; ++ipar) // loop over npar-ninitial particles from rambo
          for (int ip4 = 0; ip4 < np4; ++ip4)
            momenta[ievt][ipar+ninitial][ip4] = p_rambo[ipar][ip4];
        // Event weight
        wgts[ievt] = wgt;
        // Match the 'new double[4]' in rambo
        for (int ipar = 0; ipar < npar-ninitial; ++ipar) // loop over npar-ninitial particles from rambo
          delete[] p_rambo[ipar];
      }
    }
    return;
  }
}

std::vector<double *>  // output is a struct: momenta[npar-ninitial][4]
rambo(double et, const std::vector<double> &xm, double &wt) {
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
    xmt = xmt + abs(xm[i]);
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
    if (abs(f0) <= accu)
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
