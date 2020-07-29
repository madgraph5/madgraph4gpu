//#undef RAMBO_USES_SOA
#define RAMBO_USES_SOA 1 // presently fails check_sa.cu "error: expression must have a constant value"

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta
// Draw random momenta and the corresponding weights for nevt events
// Both initial-state and final-state particle momenta and masses are considered
void get_momenta( const int ninitial,    // input: #particles_initial
                  const double energy,   // input: energy
                  const double masses[], // input: masses[npar]
#ifdef RAMBO_USES_SOA
                  double momenta[],      // output: momenta[npar][4][nevt] as a SOA
#else
                  double momenta[],      // output: momenta[nevt][npar][4] as an AOS
#endif
                  double wgts[],         // output: weights[nevt]
                  const int npar,        // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt );      // input: #events

// Draw random momenta and the corresponding weight for a single event
// Only final-state particle momenta and masses are considered
void rambo( const double energy,    // input: energy
            const double massesf[], // input: masses[nparf]
            double p[][4],          // output: momenta[nparf][4] as struct
            double& wt,             // output: weight
            const int nparf );      // input: #particles_final (==nfinal==nexternal-ninitial)
