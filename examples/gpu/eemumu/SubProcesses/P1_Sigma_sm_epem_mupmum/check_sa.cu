#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

#include "CPPProcess.h"
// sr fixme // because of this include this needs to be a cuda file...

#include "rambo.h"

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return strlen(s) == t - s;
}

int usage(int ret = 0) {
  std::cout << "call me correctly" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false, debug = false;
  int numevts = 0;
  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0 || strcmp(argv[argn], "-v") == 0)
      verbose = true;
    else if (strcmp(argv[argn], "--debug") == 0 ||
             strcmp(argv[argn], "-d") == 0)
      debug = true;
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
  CPPProcess process(verbose, debug);

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  double energy = 1500;
  double weight;

  for (int x = 0; x < numevts; ++x) {

    // Get phase space point
    std::vector<double *> p =
        get_momenta(process.ninitial, energy, process.getMasses(), weight);

    // Set momenta for this event
    process.setMomenta(p);

    // Evaluate matrix element
    process.sigmaKin();

    const double *matrix_elements = process.getMatrixElements();

    if (verbose) {
      std::cout << "Momenta:" << std::endl;
      for (int i = 0; i < process.nexternal; i++)
        std::cout << std::setw(4) << i + 1 << setiosflags(std::ios::scientific)
                  << std::setw(14) << p[i][0]
                  << setiosflags(std::ios::scientific) << std::setw(14)
                  << p[i][1] << setiosflags(std::ios::scientific)
                  << std::setw(14) << p[i][2]
                  << setiosflags(std::ios::scientific) << std::setw(14)
                  << p[i][3] << std::endl;
      std::cout << "-----------------------------------------------------------"
                   "------------------"
                << std::endl;
    }

    if (verbose) {
      // Display matrix elements
      for (int i = 0; i < process.nprocesses; i++)
        std::cout << " Matrix element = "
                  //	 << setiosflags(ios::fixed) << setprecision(17)
                  << matrix_elements[i] << " GeV^"
                  << -(2 * process.nexternal - 8) << std::endl;

      std::cout << "-----------------------------------------------------------"
                   "------------------"
                << std::endl;
    }
  }
}
