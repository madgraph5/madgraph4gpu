#include <vector>

std::vector<std::vector<double *>> 
get_momenta(int ninitial, double energy, const std::vector<double> masses, double &wgt, int nevt);

std::vector<double *> 
rambo(double et, const std::vector<double> &xm, double &wt);
