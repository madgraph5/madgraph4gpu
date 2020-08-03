#include "mgOnGpuConfig.h"

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta
// Draw random momenta and the corresponding weights for nevt events
// Both initial-state and final-state particle momenta and masses are considered
// The number of final-state particles is nparf = npar - ninitial
void get_momenta( const int ninitial,     // input: #particles_initial
                  const double energy,    // input: energy
                  const double masses[],  // input: masses[npar]
#if defined MGONGPU_LAYOUT_ASA
                  const double rnarray[], // input: randomnumbers in [0,1] as AOSOA[npag][nparf][4][nepp]
                  double momenta1d[],     // output: momenta as AOSOA[npag][npar][4][nepp] 
#elif defined MGONGPU_LAYOUT_SOA
                  const double rnarray[], // input: randomnumbers in [0,1] as SOA[nparf][4][nevt]
                  double momenta1d[],     // output: momenta as SOA[npar][4][nevt] 
#elif defined MGONGPU_LAYOUT_AOS
                  const double rnarray[], // input: randomnumbers in [0,1] as SOA[nevt][nparf][4]
                  double momenta1d[],     // output: momenta as AOS[nevt][npar][4]
#endif
                  double wgts[],          // output: weights[nevt]
                  const int npar,         // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt );       // input: #events (==ndim=gputhread*gpublocks)

// Draw random momenta and the corresponding weight for nevt events.
// *** NB: vrambo only uses final-state masses and fills in final-state momenta,
// *** however the input masses array and output momenta array include initial-state particles
// The number of final-state particles is nparf = npar - ninitial
void vrambo( const int ninitial,       // input: #particles_initial
             const double energy,      // input: energy
             const double masses[],    // input: masses[npar] 
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
             const int npar,           // input: #particles (==nexternal==nfinal+ninitial)
             const int nevt );         // input: #events

// Generate the random numbers needed to process nevt events in rambo
#if defined MGONGPU_LAYOUT_ASA
// AOSOA: rnarray[npag][nparf][np4][nepp] where nevt=npag*nepp
#elif defined MGONGPU_LAYOUT_SOA
// SOA: rnarray[nparf][np4][nevt]
#elif defined MGONGPU_LAYOUT_AOS
// AOS: rnarray[nevt][nparf][np4]
#endif
void generateRnArray( double rnarray[], // output: randomnumbers in [0,1]
                      const int nparf,  // input: #particles_final
                      const int nevt ); // input: #events
