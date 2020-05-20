#include <algorithm>
#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>

#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include "rambo.h"
#include "timer.h"

#define gpuErrchk3(ans)                                                        \
  { gpuAssert3((ans), __FILE__, __LINE__); }

inline void gpuAssert3(cudaError_t code, const char *file, int line,
                       bool abort = true) {
  if (code != cudaSuccess) {
    printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
  }
}

#define TIMERTYPE std::chrono::high_resolution_clock

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
  bool verbose = false, debug = false, perf = false;
  int numiter = 0, gpublocks = 1, gputhreads = 1;
  std::vector<int> numvec;
  Timer<TIMERTYPE> timer;
  std::vector<float> wavetimes;

  // remove
  int nprocesses = 1;

  for (int argn = 1; argn < argc; ++argn) {
    if (strcmp(argv[argn], "--verbose") == 0 || strcmp(argv[argn], "-v") == 0)
      verbose = true;
    else if (strcmp(argv[argn], "--debug") == 0 ||
             strcmp(argv[argn], "-d") == 0)
      debug = true;
    else if (strcmp(argv[argn], "--performance") == 0 ||
             strcmp(argv[argn], "-p") == 0)
      perf = true;
    else if (is_number(argv[argn]))
      numvec.push_back(atoi(argv[argn]));
    else
      return usage(1);
  }
  int veclen = numvec.size();
  if (veclen == 3) {
    gpublocks = numvec[0];
    gputhreads = numvec[1];
    numiter = numvec[2];
  } else if (veclen == 1) {
    numiter = numvec[0];
  } else {
    return usage(1);
  }

  if (numiter == 0)
    return usage(1);

  if (verbose)
    std::cout << "# iterations: " << numiter << std::endl;

  // Create a process object
  CPPProcess process(numiter, gpublocks, gputhreads, verbose, debug);

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  double energy = 1500;
  double weight;

  int dim = gpublocks * gputhreads;

  // GPU memory
  double(*tp)[4][4];
  gpuErrchk3(cudaMalloc(&tp, dim * 4 * 4 * sizeof(double)));

  double(*matrix_element)[1];
  gpuErrchk3(
      cudaMallocManaged(&matrix_element, dim * nprocesses * sizeof(double)));

  std::vector<double> matrixelementvector;

  // later
  /*
    // Local variables and constants
    const int ncomb = 16;
    static bool goodhel[ncomb] = {ncomb * false};
    static int ntry = 0, sum_hel = 0, ngood = 0;
    static int igood[ncomb];
    static int jhel;
  */

  for (int x = 0; x < numiter; ++x) {

    // Get phase space point
    std::vector<std::vector<double *>> p =
        get_momenta(process.ninitial, energy, process.getMasses(), weight, dim);

    // Set momenta for this event
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < 4; ++i) {
        gpuErrchk3(cudaMemcpy((void *)tp[d][i], (void *)p[d][i],
                              4 * sizeof(double), cudaMemcpyHostToDevice));
      }
    }
    process.preSigmaKin();

    if (perf) {
      timer.Start();
    }

    // Evaluate matrix element
    // later process.sigmaKin(ncomb, goodhel, ntry, sum_hel, ngood, igood,
    // jhel);
    gMG5_sm::sigmaKin<<<gpublocks, gputhreads>>>(tp, matrix_element, debug,
                                                 verbose);
    cudaDeviceSynchronize();

    if (perf) {
      float gputime = timer.GetDuration();
      wavetimes.push_back(gputime);
      if (verbose)
        std::cout << "Wave function time: " << gputime << std::endl;
    }

    if (verbose || perf) {

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
          std::cout << std::string(80, '-') << std::endl;
        }

        // Display matrix elements
        for (int i = 0; i < process.nprocesses; i++) {
          if (verbose)
            std::cout << " Matrix element = "
                      //	 << setiosflags(ios::fixed) << setprecision(17)
                      << matrix_element[d][i] << " GeV^"
                      << -(2 * process.nexternal - 8) << std::endl;
          if (perf)
            matrixelementvector.push_back(matrix_element[d][i]);
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }
    } else if (!debug) {
      std::cout << ".";
    }

    for (std::vector<std::vector<double *>>::iterator it = p.begin();
         it != p.end(); ++it) {
      for (std::vector<double *>::iterator jt = it->begin(); jt != it->end();
           ++jt) {
        delete[] & (**jt);
      }
    }
  }

  if (!(verbose || debug || perf)) {
    std::cout << std::endl;
  }

  if (perf) {

    float sum = std::accumulate(wavetimes.begin(), wavetimes.end(), 0.0);
    int num_wts = wavetimes.size();
    float mean = sum / num_wts;
    float sq_sum = std::inner_product(wavetimes.begin(), wavetimes.end(),
                                      wavetimes.begin(), 0.0);
    float stdev = std::sqrt(sq_sum / num_wts - mean * mean);
    std::vector<float>::iterator mintime =
        std::min_element(wavetimes.begin(), wavetimes.end());
    std::vector<float>::iterator maxtime =
        std::max_element(wavetimes.begin(), wavetimes.end());

    int num_mes = matrixelementvector.size();
    std::vector<double>::iterator maxelem = std::max_element(
        matrixelementvector.begin(), matrixelementvector.end());
    std::vector<double>::iterator minelem = std::min_element(
        matrixelementvector.begin(), matrixelementvector.end());

    std::cout << "***********************************" << std::endl
              << "NumIterations        = " << numiter << std::endl
              << "NumThreadsPerBlock   = " << gputhreads << std::endl
              << "NumBlocksPerGrid     = " << gpublocks << std::endl
              << "NumberOfEntries      = " << num_wts << std::endl
              << std::scientific << "TotalTimeInWaveFuncs = " << sum
              << std::endl
              << "MeanTimeinWaveFuncs  = " << mean << std::endl
              << "StdDevWaveFuncs      = " << stdev << std::endl
              << "MinTimeInWaveFuncs   = " << *mintime << std::endl
              << "MaxTimeInWaveFuncs   = " << *maxtime << std::endl
              << "-----------------------------------" << std::endl
              << "ProcessID:           = " << getpid() << std::endl
              << "NProcesses           = " << process.nprocesses << std::endl
              << "NumMatrixElements    = " << num_mes << std::endl
              << std::scientific << "MaxMatrixElemValue   = " << *maxelem
              << std::endl
              << "MinMatrixElemValue   = " << *minelem << std::endl;
  }
}
