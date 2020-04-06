#include <cstring>
#include <iomanip>
#include <iostream>

#include "CPPProcess.h"
#include "rambo.h"

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return strlen(s) == t - s;
}

int usage(int ret = 0) {
  std::cout << "do it correctly" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false;
  int numevts = 0;
  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0)
      verbose = true;
    else if (is_number(argv[argn]))
      numevts = atoi(argv[argn]);
    else
      return usage(1);
  }
  if (numevts == 0)
    return usage(1);

  if (verbose)
    std::cout << "num evts: " << numevts << std::endl;

  // Create a process object
  CPPProcess process;

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat", verbose);

  double energy = 1500;
  double weight;

  for (int x = 0; x <= numevts; ++x) {

    // Get phase space point
    vector<double *> p =
        get_momenta(process.ninitial, energy, process.getMasses(), weight);

    // Set momenta for this event
    process.setMomenta(p);

    // Evaluate matrix element
    process.sigmaKin(verbose);

    const double *matrix_elements = process.getMatrixElements();

    if (verbose) {
      cout << "Momenta:" << endl;
      for (int i = 0; i < process.nexternal; i++)
        cout << setw(4) << i + 1 << setiosflags(ios::scientific) << setw(14)
             << p[i][0] << setiosflags(ios::scientific) << setw(14) << p[i][1]
             << setiosflags(ios::scientific) << setw(14) << p[i][2]
             << setiosflags(ios::scientific) << setw(14) << p[i][3] << endl;
      cout << " ---------------------------------------------------------------"
              "--------------"
           << endl;
    }

    if (verbose) {
      // Display matrix elements
      for (int i = 0; i < process.nprocesses; i++)
        cout << " Matrix element = "
             //	 << setiosflags(ios::fixed) << setprecision(17)
             << matrix_elements[i] << " GeV^" << -(2 * process.nexternal - 8)
             << endl;

      cout << " ---------------------------------------------------------------"
              "--------------"
           << endl;
    }
  }

  process.printPerformanceStats();
}
