#ifndef RAMBO2TONM0_H
#define RAMBO2TONM0_H 1

#include "mgOnGpuConfig.h"

// Simplified rambo version for 2 to N (with N>=2) processes with massless particles
namespace rambo2toNm0
{
  /*
  const int np4 = 4; // the dimension of 4-momenta (E,px,py,pz)

  const int npari = 2; // #particles in the initial state
  const int nparf = 2; // #particles in the final state
  const int npar = npari + nparf; // #particles in total
  */

  // Fill in the momenta of the initial particles
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  void getMomentaInitial( const double energy,    // input: energy
#if defined MGONGPU_LAYOUT_ASA
                          double momenta1d[],     // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                          double momenta1d[],     // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                          double momenta1d[],     // output: momenta as AOS[nevt][npar][4]
#endif
                          const int nevt );       // input: #events

  // Fill in the momenta of the final particles using the RAMBO algorithm
  // [NB: the output buffer includes both initial and final momenta, but only initial momenta are filled in]
  void getMomentaFinal( const double energy,    // input: energy
#if defined MGONGPU_LAYOUT_ASA
                        const double rnarray1d[], // input: randomnumbers in [0,1] as AOSOA[npag][nparf][4][nepp]
                        double momenta1d[],       // output: momenta as AOSOA[npag][npar][4][nepp]
#elif defined MGONGPU_LAYOUT_SOA
                        const double rnarray1d[], // input: randomnumbers in [0,1] as SOA[nparf][4][nevt]
                        double momenta1d[],       // output: momenta as SOA[npar][4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
                        const double rnarray1d[], // input: randomnumbers in [0,1] as SOA[nevt][nparf][4]
                        double momenta1d[],       // output: momenta as AOS[nevt][npar][4]
#endif
                        double wgts[],            // output: weights[nevt]
                        const int nevt );         // input: #events

  // Generate the random numbers needed to process nevt events in rambo
#if defined MGONGPU_LAYOUT_ASA
  // AOSOA: rnarray[npag][nparf][np4][nepp] where nevt=npag*nepp
#elif defined MGONGPU_LAYOUT_SOA
  // SOA: rnarray[nparf][np4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
  // AOS: rnarray[nevt][nparf][np4]
#endif
  void generateRnArray( double rnarray1d[], // output: randomnumbers in [0,1]
                        const int nevt );   // input: #events

}

#endif // RAMBO2TONM0_H 1
