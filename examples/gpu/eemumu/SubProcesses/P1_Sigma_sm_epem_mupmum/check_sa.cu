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
  int numevts = 0, gpuwarps = 1, gputhreads = 1;
  std::vector<int> numvec;
  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0 || strcmp(argv[argn], "-v") == 0)
      verbose = true;
    else if (strcmp(argv[argn], "--debug") == 0 ||
             strcmp(argv[argn], "-d") == 0)
      debug = true;
    else if (is_number(argv[argn]))
      numvec.push_back(atoi(argv[argn]));
    // numevts = atoi(argv[argn]);
    else
      return usage(1);
  }
  int veclen = numvec.size();
  if (veclen == 3) {
    gpuwarps = numvec[0];
    gputhreads = numvec[1];
    numevts = numvec[2];
  } else if (veclen == 1) {
    numevts = numvec[0];
  } else {
    return usage(1);
  }

  if (numevts == 0)
    return usage(1);

  if (verbose)
    std::cout << "num evts: " << numevts << std::endl;

  // Create a process object
  CPPProcess process(gpuwarps, gputhreads, verbose, debug);

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  double energy = 1500;
  double weight;

  int dim = process.getDim();

  for (int x = 0; x < numevts; ++x) {

    if (!(verbose || debug)) {
      std::cout << ".";
    }

    // Get phase space point
    std::vector<std::vector<double *>> p =
        get_momenta(process.ninitial, energy, process.getMasses(), weight, dim);

    // Set momenta for this event
    process.setMomenta(p);

    // Evaluate matrix element
    process.sigmaKin();

    double **matrix_elements = process.getMatrixElements();

    for (int d = 0; d < dim; ++d) {

      if (verbose) {
        std::cout << "Momenta:" << std::endl;
        for (int i = 0; i < process.nexternal; i++)
          std::cout << std::setw(4) << i + 1
                    << setiosflags(std::ios::scientific) << std::setw(14)
                    << p[d][i][0] << setiosflags(std::ios::scientific)
                    << std::setw(14) << p[d][i][1]
                    << setiosflags(std::ios::scientific) << std::setw(14)
                    << p[d][i][2] << setiosflags(std::ios::scientific)
                    << std::setw(14) << p[d][i][3] << std::endl;
        std::cout
            << "-----------------------------------------------------------"
               "------------------"
            << std::endl;
      }

      if (verbose) {
        // Display matrix elements
        for (int i = 0; i < process.nprocesses; i++)
          std::cout << " Matrix element = "
                    //	 << setiosflags(ios::fixed) << setprecision(17)
                    << matrix_elements[d][i] << " GeV^"
                    << -(2 * process.nexternal - 8) << std::endl;

        std::cout
            << "-----------------------------------------------------------"
               "------------------"
            << std::endl;
      }
    }
    for (std::vector<std::vector<double *>>::iterator it = p.begin();
         it != p.end(); ++it) {
      for (std::vector<double *>::iterator jt = it->begin(); jt != it->end();
           ++jt) {
        delete[] & (**jt);
      }
    }
  }

  if (!(verbose || debug)) {
    std::cout << std::endl;
  }
}
