//#include <vector>

//std::vector<std::vector<double *>> // output is an AOS: momenta[nevt][nexternal][4]
//get_momenta(int ninitial, double energy, const std::vector<double> masses, double &wgt, int nevt);

// Auxiliary function to change convention between MadGraph5_aMC@NLO and rambo four momenta
// Draw random momenta and the corresponding weights for nevt events
// Both initial-state and final-state particle momenta and masses are considered
void get_momenta( const int ninitial,    // input: #particles_initial
                  const double energy,   // input: energy
                  const double masses[], // input: masses[npar]
                  double momenta[],      // output: momenta[nevt][npar][4] as an AOS
                  double wgts[],         // output: weights[nevt]
                  const int npar,        // input: #particles (==nexternal==nfinal+ninitial)
                  const int nevt );      // input: #events

//std::vector<double *>  // output is a struct: momenta[npar-ninitial][4]
//rambo(double et, const std::vector<double> &xm, double &wt);

// Draw random momenta and the corresponding weight for a single event
// Only final-state particle momenta and masses are considered
void rambo( const double energy,    // input: energy
            const double massesf[], // input: masses[nparf]
            double p[][4],          // output: momenta[nparf][4] as struct
            double& wt,             // output: weight
            const int nparf );      // input: #particles_final (==nfinal==nexternal-ninitial)
