#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

#include "CPPProcess.h"
// sr fixme // because of this include this needs to be a cuda file...

#include "rambo.h"

// see https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
template <typename T>
T **alloc2DArray(int y, int x, int &size = 0, bool init = false) {

  int i, j, count = 0;
  T *ptr, **arr;

  size = sizeof(T *) * y + sizeof(T) * y * x;
  arr = (T **)malloc(size);

  ptr = (T *)(arr + y); // points to first T value

  // for loop to point rows pointer to appropriate location in 2D array
  for (i = 0; i < y; ++i)
    arr[i] = (ptr + x * i);

  if (init)
    for (i = 0; i < y; i++)
      for (j = 0; j < x; j++)
        arr[i][j] = ++count;

  return arr;
}

// inspiration from 4) in
// https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
template <typename T>
T ***alloc3DArray(int z, int y, int x, int &size, bool init = false) {

  int i, j, k;
  T *ptr2t, **ptr2p, ***arr;

  size = (z + z * y) * sizeof(T *) + z * y * x * sizeof(T);
  arr = (T ***)malloc(size);

  ptr2p = (T **)(arr + z);    // points to first 2d array
  ptr2t = (T *)(arr + z + y); // points to first T value

  // for loop to point dimension ids to array locations (3d->2d and 2d->T)
  for (i = 0; i < z; ++i) {
    arr[i] = (T **)((T *)(ptr2p + y * i) + y * x * i);
    for (j = 0; j < y; ++j) {
      arr[i][j] = ptr2t + x * j;
    }
    ptr2t = (T *)((T **)ptr2t + y) + y * x;
  }

  if (init)
    for (i = 0; i < z; ++i)
      for (j = 0; j < y; ++j)
        for (k = 0; k < x; ++k)
          arr[i][j][k] = 0;

  return arr;
}

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

  double energy = 1500, weight, ***p;
  int psize = 0;
  int dim = process.getDim();
  int nioparticles = process.getNIOParticles();

  for (int x = 0; x < numevts; ++x) {

    if (!(verbose || debug)) {
      std::cout << ".";
    }

    p = alloc3DArray<double>(
        dim, 4, 4, psize); // sr fixme // first 4 paramater should be a variable
    // Get phase space point
    get_momenta(process.ninitial, energy, process.getMasses(), weight, dim, p);

    // Set momenta for this event
    process.setMomenta(p, psize);

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
  }
  if (!(verbose || debug)) {
    std::cout << std::endl;
  }
}
