#include "mgOnGpuConfig.h" // check if RAMBO_USES_AOS

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta
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
                  double wgts[],          // output: weights[nevt]
                  const int npar,         // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt );       // input: #events

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
             const int ievt );         // input: event ID to be written out out of #events

// Generate the random numbers needed to process nevt events in rambo
#ifndef RAMBO_USES_AOS
void generateRnArray( double rnarray[], // output: randomnumbers[nparf][4][nevt] in [0,1] as an SOA
                      const int nparf,  // input: #particles_final
                      const int nevt ); // input: #events
#else
void generateRnArray( double rnarray[], // output: randomnumbers[nevt][nparf][4] in [0,1] as an AOS
                      const int nparf,  // input: #particles_final
                      const int nevt ); // input: #events
#endif
