#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "vrambo.h"
#include "Random.h"

//std::vector<std::vector<double *>> // output is an AOS: momenta[nevt][nexternal][4]
//get_momenta(int ninitial, double energy, const std::vector<double> masses, double &wgt, int nevt) {

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta.
// Draw random momenta and the corresponding weights for nevt events
// Both initial-state and final-state particle momenta and masses are considered
void get_momenta( const int ninitial,    // input: #particles_initial
                  const double energy,   // input: energy
                  const double masses[], // input: masses[npar]
                  double momenta1d[],    // output: momenta[nevt][npar][4] as an AOS
                  double wgts[],         // output: wgts[nevt]
                  const int npar,        // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt )       // input: #events
{
  const int nparf = npar - ninitial; // (previously was nfinal = nexternal - ninitial)
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
      const double* massesf = masses+ninitial; // skip the first ninitial masses
      double pf_rambo[nparf][np4]; // rambo draws random momenta only for final-state particles
      double wgt;
      rambo( m1, massesf, pf_rambo, wgt, nparf ); // NB input 'energy' is ignored for ninitial==1
      for (int iparf = 0; iparf < nparf; ++iparf) // loop over npar-ninitial particles from rambo
        for (int ip4 = 0; ip4 < np4; ++ip4)
          momenta[ievt][iparf+ninitial][ip4] = pf_rambo[iparf][ip4];
      // Event weight
      wgts[ievt] = wgt;
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
    const double m2 = masses[1];
    const double mom =
      sqrt((pow(e2, 2) - 2 * e2 * pow(m1, 2) + pow(m1, 4) -
            2 * e2 * pow(m2, 2) - 2 * pow(m1, 2) * pow(m2, 2) + pow(m2, 4)) /
           (4 * e2));
    const double energy1 = sqrt(pow(mom, 2) + pow(m1, 2));
    const double energy2 = sqrt(pow(mom, 2) + pow(m2, 2));
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
      if (nparf == 1) {
        // Momenta for the outgoing particle
        momenta[ievt][2][0] = m1;
        momenta[ievt][2][1] = 0;
        momenta[ievt][2][2] = 0;
        momenta[ievt][2][3] = 0;
        // Event weight
        wgts[ievt] = 1;
      }
      // #Initial==2, #Final>1
      else {
        // Momenta for the outgoing particles
        const double* massesf = masses+ninitial; // skip the first ninitial masses
        double pf_rambo[nparf][np4]; // rambo draws random momenta only for final-state particles
        double wgt;
        rambo( energy, massesf, pf_rambo, wgt, nparf );
        for (int iparf = 0; iparf < nparf; ++iparf) // loop over npar-ninitial particles from rambo
          for (int ip4 = 0; ip4 < np4; ++ip4)
            momenta[ievt][iparf+ninitial][ip4] = pf_rambo[iparf][ip4];
        // Event weight
        wgts[ievt] = wgt;
      }
    }
    return;
  }
}

//std::vector<double *>  // output is a struct: momenta[npar-ninitial][4]
//rambo(double et, const std::vector<double> &xm, double &wt) {

// Draw random momenta and the corresponding weight for a single event
// Only final-state particle momenta and masses are considered
void rambo( const double energy, // input: energy
            const double xmf[],  // input: masses[nparf]
            double p[][4],       // output: momenta[nparf][4] as struct
            double& wt,          // output: weight
            const int nparf )    // input: #particles_final (==nfinal==nexternal-ninitial)
{
  /****************************************************************************
   *                       rambo                                              *
   *    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)                       *
   *                                                                          *
   *    a democratic multi-particle phase space generator                     *
   *    authors:  s.d. ellis,  r. kleiss,  w.j. stirling                      *
   *    this is version 1.0 -  written by r. kleiss                           *
   *    -- adjusted by hans kuijf, weights are logarithmic (1990-08-20)       *
   *    -- adjusted by andrea valassi, remove std::vectors (2020-07-29)       *
   *                                                                          *
   *    nparf  = number of final-state particles ( =nexternal-nincoming )     *
   *    energy = total centre-of-mass energy                                  *
   *    xmf    = final-state particle masses ( dim=nexternal-nincoming )      *
   *    p      = final-state particle momenta ( dim=(nexternal-nincoming,4) ) *
   *    wt     = weight of the event                                          *
   ****************************************************************************/  
  double q[nparf][4];
  double z[nparf], r[4], b[3], p2[nparf], xmf2[nparf], e[nparf], v[nparf];
  static std::vector<int> iwarn(5, 0);
  static double acc = 1e-14;
  static int itmax = 6, ibegin = 0;
  static double twopi = 8. * atan(1.);
  static double po2log = log(twopi / 4.);

  // initialization step: factorials for the phase space weight
  if (ibegin == 0) {
    ibegin = 1;
    z[1] = po2log;
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = z[kpar - 1] + po2log - 2. * log(double(kpar - 1));
    for (int kpar = 2; kpar < nparf; kpar++)
      z[kpar] = (z[kpar] - log(double(kpar)));
  }
  // check on the number of particles
  if (nparf < 1 || nparf > 101) {
    std::cout << "Too few or many particles: " << nparf << std::endl;
    exit(-1);
  }
  // check whether total energy is sufficient; count nonzero masses
  double xmft = 0.;
  int nm = 0;
  for (int iparf = 0; iparf < nparf; iparf++) {
    if (xmf[iparf] != 0.)
      nm = nm + 1;
    xmft = xmft + abs(xmf[iparf]);
  }
  if (xmft > energy) {
    std::cout << "Too low energy: " << energy << " needed " << xmft << std::endl;
    exit(-1);
  }
  // the parameter values are now accepted

  // generate n massless momenta in infinite phase space
  for (int iparf = 0; iparf < nparf; iparf++) {
    double r1 = rn(1);
    double c = 2. * r1 - 1.;
    double s = sqrt(1. - c * c);
    double f = twopi * rn(2);
    r1 = rn(3);
    double r2 = rn(4);
    q[iparf][0] = -log(r1 * r2);
    q[iparf][3] = q[iparf][0] * c;
    q[iparf][2] = q[iparf][0] * s * cos(f);
    q[iparf][1] = q[iparf][0] * s * sin(f);
  }
  // calculate the parameters of the conformal transformation
  for (int i4 = 0; i4 < 4; i4++)
    r[i4] = 0.;
  for (int iparf = 0; iparf < nparf; iparf++) {
    for (int i4 = 0; i4 < 4; i4++)
      r[i4] = r[i4] + q[iparf][i4];
  }
  double rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
  for (int i4 = 1; i4 < 4; i4++)
    b[i4-1] = -r[i4] / rmas;
  double g = r[0] / rmas;
  double a = 1. / (1. + g);
  double x = energy / rmas;

  // transform the q's conformally into the p's
  for (int iparf = 0; iparf < nparf; iparf++) {
    double bq = b[0] * q[iparf][1] + b[1] * q[iparf][2] + b[2] * q[iparf][3];
    for (int i4 = 1; i4 < 4; i4++)
      p[iparf][i4] = x * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
    p[iparf][0] = x * (g * q[iparf][0] + bq);
  }

  // calculate weight and possible warnings (NB return log of weight)
  wt = po2log;
  if (nparf != 2)
    wt = (2. * nparf - 4.) * log(energy) + z[nparf-1];
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
    return;
  }

  // massive particles: rescale the momenta by a factor x
  double xmax = sqrt(1. - pow(xmft / energy, 2));
  for (int iparf = 0; iparf < nparf; iparf++) {
    xmf2[iparf] = pow(xmf[iparf], 2);
    p2[iparf] = pow(p[iparf][0], 2);
  }
  int iter = 0;
  x = xmax;
  double accu = energy * acc;
  while (true) {
    double f0 = -energy;
    double g0 = 0.;
    double x2 = x * x;
    for (int iparf = 0; iparf < nparf; iparf++) {
      e[iparf] = sqrt(xmf2[iparf] + x2 * p2[iparf]);
      f0 = f0 + e[iparf];
      g0 = g0 + p2[iparf] / e[iparf];
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
  for (int iparf = 0; iparf < nparf; iparf++) {
    v[iparf] = x * p[iparf][0];
    for (int i4 = 1; i4 < 4; i4++)
      p[iparf][i4] = x * p[iparf][i4];
    p[iparf][0] = e[iparf];
  }

  // calculate the mass-effect weight factor
  double wt2 = 1.;
  double wt3 = 0.;
  for (int iparf = 0; iparf < nparf; iparf++) {
    wt2 = wt2 * v[iparf] / e[iparf];
    wt3 = wt3 + pow(v[iparf], 2) / e[iparf];
  }
  double wtm = (2. * nparf - 3.) * log(x) + log(wt2 / wt3 * energy);

  // return for weighted massive momenta
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
  return;
}
