#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

#include "CPPProcess.h"
// sr fixme // because of this include this needs to be a cuda file...

#include "rambo.h"
//#include <thrust/host_vector.h>

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

  processMem *pm = new processMem();

  cudaDeviceSynchronize();

  cudaMallocManaged(&pm->tamp,
                    pm->tnamplitudes * sizeof(thrust::complex<double>));

  cudaMallocManaged(&pm->tw, pm->twrows * sizeof(thrust::complex<double>));
  // cudaDeviceSynchronize();
  for (int i = 0; i < pm->twrows; ++i) {
    cudaMallocManaged(&pm->tw[i],
                      pm->tnwavefuncs * sizeof(thrust::complex<double>));
    // cudaDeviceSynchronize();
    for (int y = 0; y < pm->twrows; ++y) {
      pm->tw[i][y] =
          thrust::complex<double>(0.00, 0.00); // double(i), double(y));
    }
  }

  cudaDeviceSynchronize();
  // Create a process object
  CPPProcess process(pm);

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat", verbose);

  double energy = 1500;
  double weight;

  for (int x = 0; x < numevts; ++x) {

    process.resetGPUMemory();

    // Get phase space point
    std::vector<double *> p =
        get_momenta(process.ninitial, energy, process.getMasses(), weight);

    // Set momenta for this event
    process.setMomenta(p);

    // Evaluate matrix element
    process.sigmaKin(verbose);

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
