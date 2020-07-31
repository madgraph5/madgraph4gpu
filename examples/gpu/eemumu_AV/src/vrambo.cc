#include <cmath>
#include <cstdlib>
#include <iostream>

#include "vrambo.h"
#include "Random.h"

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta.
// Draw random momenta and the corresponding weights for nevt events
// Both initial-state and final-state particle momenta and masses are considered
// The number of final-state particles is nparf = npar - ninitial
void get_momenta( const int ninitial,     // input: #particles_initial
                  const double energy,    // input: energy
                  const double masses[],  // input: masses[npar]
#ifndef RAMBO_USES_AOS
                  const double rnarray[], // input: randomnumbers[nparf][4][nevt] in [0,1] as an SOA
                  double momenta1d[],     // output: momenta[npar][4][nevt] as a SOA
#else
                  const double rnarray[], // input: randomnumbers[nevt][nparf][4] in [0,1] as an AOS
                  double momenta1d[],     // output: momenta[nevt][npar][4] as an AOS
#endif
                  double wgts[],          // output: wgts[nevt]
                  const int npar,         // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt )        // input: #events
{
  const int nparf = npar - ninitial; // (previously was nfinal = nexternal - ninitial)
  const double e2 = pow(energy, 2);
  const double m1 = masses[0];

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
#ifndef RAMBO_USES_AOS
  double (*momenta)[np4][nevt] = (double (*)[np4][nevt]) momenta1d; // cast to multiD array pointer (SOA)
#else
  double (*momenta)[npar][np4] = (double (*)[npar][np4]) momenta1d; // cast to multiD array pointer (AOS)
#endif

  // #Initial==1
  if (ninitial == 1) {
    for (int ievt = 0; ievt < nevt; ++ievt) {
      // Momenta for the incoming particle
#ifndef RAMBO_USES_AOS
      momenta[0][0][ievt] = m1;
      momenta[0][1][ievt] = 0;
      momenta[0][2][ievt] = 0;
      momenta[0][3][ievt] = 0;
#else
      momenta[ievt][0][0] = m1;
      momenta[ievt][0][1] = 0;
      momenta[ievt][0][2] = 0;
      momenta[ievt][0][3] = 0;
#endif
      // Momenta for the outgoing particles and event weight
      vrambo( ninitial, m1, // NB input 'energy' is ignored for ninitial==1
              masses, rnarray, (double*)momenta, wgts, npar, nevt, ievt );
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
#ifndef RAMBO_USES_AOS
      momenta[0][0][ievt] = energy1;
      momenta[0][1][ievt] = 0;
      momenta[0][2][ievt] = 0;
      momenta[0][3][ievt] = mom;
      momenta[1][0][ievt] = energy2;
      momenta[1][1][ievt] = 0;
      momenta[1][2][ievt] = 0;
      momenta[1][3][ievt] = -mom;
#else
      momenta[ievt][0][0] = energy1;
      momenta[ievt][0][1] = 0;
      momenta[ievt][0][2] = 0;
      momenta[ievt][0][3] = mom;
      momenta[ievt][1][0] = energy2;
      momenta[ievt][1][1] = 0;
      momenta[ievt][1][2] = 0;
      momenta[ievt][1][3] = -mom;
#endif
      // #Initial==2, #Final==1
      if (nparf == 1) {
        // Momenta for the outgoing particle
#ifndef RAMBO_USES_AOS
        momenta[2][0][ievt] = m1;
        momenta[2][1][ievt] = 0;
        momenta[2][2][ievt] = 0;
        momenta[2][3][ievt] = 0;
#else
        momenta[ievt][2][0] = m1;
        momenta[ievt][2][1] = 0;
        momenta[ievt][2][2] = 0;
        momenta[ievt][2][3] = 0;
#endif
        // Event weight
        wgts[ievt] = 1;
      }
      // #Initial==2, #Final>1
      else {
        // Momenta for the outgoing particles and event weight
        vrambo( ninitial, energy, masses, rnarray, (double*)momenta, wgts, npar, nevt, ievt );
      }
    }
    return;
  }
}

// Draw random momenta and the corresponding weight for event ievt out of nevt
// *** NB: vrambo only uses final-state masses and fills in final-state momenta,
// *** however the input masses array and output momenta array include initial-state particles
// The number of final-state particles is nparf = npar - ninitial
void vrambo( const int ninitial,       // input: #particles_initial
             const double energy,      // input: energy
             const double masses[],    // input: masses[npar] 
#ifndef RAMBO_USES_AOS
             const double rnarray1d[], // input: randomnumbers[nparf][4][nevt] in [0,1] as an SOA
             double momenta1d[],       // output: momenta[npar][4][nevt] as a SOA
#else
             const double rnarray1d[], // input: randomnumbers[nevt][nparf][4] in [0,1] as an AOS
             double momenta1d[],       // output: momenta[nevt][npar][4] as an AOS
#endif
             double wgts[],            // output: weights[nevt]
             const int npar,           // input: #particles (==nexternal==nfinal+ninitial)
             const int nevt,           // input: #events
             const int ievt )          // input: event ID to be written out out of #events
{
  /****************************************************************************
   *                       rambo                                              *
   *    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)                       *
   *                                                                          *
   *    a democratic multi-particle phase space generator                     *
   *    authors:  s.d. ellis,  r. kleiss,  w.j. stirling                      *
   *    this is version 1.0 -  written by r. kleiss                           *
   *    -- adjusted by hans kuijf, weights are logarithmic (1990-08-20)       *
   *    -- adjusted by madgraph@sheffield_gpu_hackathon team (2020-07-29)     *
   *                                                                          *
   *    nparf  = number of final-state particles ( =nexternal-nincoming )     *
   *    energy = total centre-of-mass energy                                  *
   *    xmf    = final-state particle masses ( dim=nexternal-nincoming )      *
   *    p      = final-state particle momenta ( dim=(nexternal-nincoming,4) ) *
   *    wt     = weight of the event                                          *
   ****************************************************************************/  
  const int nparf = npar - ninitial;

  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)
#ifndef RAMBO_USES_AOS
  double (*rnarray)[np4][nevt] = (double (*)[np4][nevt]) rnarray1d; // cast to multiD array pointer (SOA)
  double (*momenta)[np4][nevt] = (double (*)[np4][nevt]) momenta1d; // cast to multiD array pointer (SOA)
#else
  double (*rnarray)[nparf][np4] = (double (*)[nparf][np4]) rnarray1d; // cast to multiD array pointer (AOS)
  double (*momenta)[npar][np4] = (double (*)[npar][np4]) momenta1d; // cast to multiD array pointer (AOS)
#endif
  double& wt = wgts[ievt];
  const double* xmf = masses+ninitial; // skip the first ninitial masses

  double q[nparf][np4];
  double z[nparf], r[np4], b[np4-1], p2[nparf], xmf2[nparf], e[nparf], v[nparf];
  int iwarn[5] = {0,0,0,0,0};
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
#ifndef RAMBO_USES_AOS
    const double r1 = rnarray[iparf][0][ievt];
    const double r2 = rnarray[iparf][1][ievt];
    const double r3 = rnarray[iparf][2][ievt];
    const double r4 = rnarray[iparf][3][ievt];
#else
    const double r1 = rnarray[ievt][iparf][0];
    const double r2 = rnarray[ievt][iparf][1];
    const double r3 = rnarray[ievt][iparf][2];
    const double r4 = rnarray[ievt][iparf][3];
#endif
    const double c = 2. * r1 - 1.;
    const double s = sqrt(1. - c * c);
    const double f = twopi * r2;
    q[iparf][0] = -log(r3 * r4);
    q[iparf][3] = q[iparf][0] * c;
    q[iparf][2] = q[iparf][0] * s * cos(f);
    q[iparf][1] = q[iparf][0] * s * sin(f);
  }
  // calculate the parameters of the conformal transformation
  for (int i4 = 0; i4 < np4; i4++)
    r[i4] = 0.;
  for (int iparf = 0; iparf < nparf; iparf++) {
    for (int i4 = 0; i4 < np4; i4++)
      r[i4] = r[i4] + q[iparf][i4];
  }
  double rmas = sqrt(pow(r[0], 2) - pow(r[3], 2) - pow(r[2], 2) - pow(r[1], 2));
  for (int i4 = 1; i4 < np4; i4++)
    b[i4-1] = -r[i4] / rmas;
  double g = r[0] / rmas;
  double a = 1. / (1. + g);
  double x = energy / rmas;

  // transform the q's conformally into the p's (i.e. the 'momenta')
  for (int iparf = 0; iparf < nparf; iparf++) {
    double bq = b[0] * q[iparf][1] + b[1] * q[iparf][2] + b[2] * q[iparf][3];
#ifndef RAMBO_USES_AOS
    for (int i4 = 1; i4 < np4; i4++)
      momenta[iparf+ninitial][i4][ievt] = x * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
    momenta[iparf+ninitial][0][ievt] = x * (g * q[iparf][0] + bq);
#else
    for (int i4 = 1; i4 < np4; i4++)
      momenta[ievt][iparf+ninitial][i4] = x * (q[iparf][i4] + b[i4-1] * (q[iparf][0] + a * bq));
    momenta[ievt][iparf+ninitial][0] = x * (g * q[iparf][0] + bq);
#endif
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
#ifndef RAMBO_USES_AOS
    p2[iparf] = pow(momenta[iparf+ninitial][0][ievt], 2);
#else
    p2[iparf] = pow(momenta[ievt][iparf+ninitial][0], 2);
#endif
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
#ifndef RAMBO_USES_AOS
    v[iparf] = x * momenta[iparf+ninitial][0][ievt];
    for (int i4 = 1; i4 < np4; i4++)
      momenta[iparf+ninitial][i4][ievt] = x * momenta[iparf+ninitial][i4][ievt];
    momenta[iparf+ninitial][0][ievt] = e[iparf];
#else
    v[iparf] = x * momenta[ievt][iparf+ninitial][0];
    for (int i4 = 1; i4 < np4; i4++)
      momenta[ievt][iparf+ninitial][i4] = x * momenta[ievt][iparf+ninitial][i4];
    momenta[ievt][iparf+ninitial][0] = e[iparf];
#endif
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

// Generate the random numbers needed to process nevt events in rambo
#ifndef RAMBO_USES_AOS
void generateRnArray( double rnarray1d[], // output: randomnumbers[nparf][4][nevt] in [0,1] as an SOA
                      const int nparf,    // input: #particles_final
                      const int nevt )    // input: #events
{
  const int np4 = 4; // 4 random numbers (like the dimension of 4-momenta) are needed for each particle
  double (*rnarray)[np4][nevt] = (double (*)[np4][nevt]) rnarray1d; // cast to multiD array pointer (AOS)
  for (int iparf = 0; iparf < nparf; iparf++)
    for (int ip4 = 0; ip4 < np4; ++ip4)
      for (int ievt = 0; ievt < nevt; ++ievt)
        rnarray[iparf][ip4][ievt] = rn(0); // SOA[nparf][np4][nevt]
}
#else
void generateRnArray( double rnarray1d[], // output: randomnumbers[nevt][nparf][4] in [0,1] as an AOS
                      const int nparf,    // input: #particles_final
                      const int nevt )    // input: #events
{
  const int np4 = 4; // 4 random numbers (like the dimension of 4-momenta) are needed for each particle
  double (*rnarray)[nparf][np4] = (double (*)[nparf][np4]) rnarray1d; // cast to multiD array pointer (SOA)
  for (int ievt = 0; ievt < nevt; ++ievt)
    for (int iparf = 0; iparf < nparf; iparf++)
      for (int ip4 = 0; ip4 < np4; ++ip4)
        rnarray[ievt][iparf][ip4] = rn(0); // AOS[nevt][nparf][np4]
}
#endif
