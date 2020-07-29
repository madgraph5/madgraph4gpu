#include <vector>

//std::vector<std::vector<double *>> // output is an AOS: momenta[nevt][nexternal][4]
//get_momenta(int ninitial, double energy, const std::vector<double> masses, double &wgt, int nevt);

void get_momenta( int ninitial,          // input: #particles_initial
                  double energy,         // input: energy
                  const double masses[], // input: masses[npar]
                  double momenta[],      // output: momenta[nevt][npar][4] as an AOS
                  double wgts[],         // output: wgts[nevt]
                  int npar,              // input: #particles (==nexternal==nfinal+ninitial)
                  int nevt );            // input: #events

std::vector<double *>  // output is a struct: momenta[npar-ninitial][4]
rambo(double et, const std::vector<double> &xm, double &wt);

//std::vector<double *>  // output is a struct: momenta[npar-ninitial][4]
//void rambo( double et, 
//            const std::vector<double> &xm, 
//            double &wt );
