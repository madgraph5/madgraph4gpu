#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>

#include "CPPProcess.h"
//#include "HelAmps_sm.h"

#include "rambo.h"
#include "timer.h"

#include "Kokkos_Core.hpp"

#define TIMERTYPE std::chrono::high_resolution_clock

bool is_number(const char *s) {
  const char *t = s;
  while (*t != '\0' && isdigit(*t))
    ++t;
  return strlen(s) == t - s;
}

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0 
            << " [--verbose|-v] [--debug|-d] [--performance|-p]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false, debug = false, perf = false;
  int numiter = 0, gpublocks = 1, gputhreads = 1;
  std::vector<int> numvec;
  Timer<TIMERTYPE> timer;
  std::vector<float> wavetimes;


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
      return usage(argv[0]);
  }
  int veclen = numvec.size();
  if (veclen == 3) {
    gpublocks = numvec[0];
    gputhreads = numvec[1];
    numiter = numvec[2];
  } else if (veclen == 1) {
    numiter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (numiter == 0)
    return usage(argv[0]);

  Kokkos::initialize(argc, argv);
  if (verbose)
    std::cout << "# iterations: " << numiter << std::endl;
  { // start Kokkos View space
    // Create a process object
    CPPProcess<Kokkos::DefaultExecutionSpace> process(numiter, gpublocks, gputhreads);

    // Read param_card and set parameters
    process.initProc("../../Cards/param_card.dat");

    double energy = 1500;
    double weight;

    int meGeVexponent = -(2 * process.nexternal - 8);
    int dim = gpublocks * gputhreads;

    Kokkos::View<double*,Kokkos::DefaultExecutionSpace> meDevPtr("meDevPtr",dim*1);
    auto meHostPtr = Kokkos::create_mirror_view(meDevPtr);

    std::vector<double> matrixelementvector;

    for (int x = 0; x < numiter; ++x) {
      // printf("iter %d of %d\n",x,numiter);
      // Get phase space point
      auto p = get_momenta(process.ninitial, process.nexternal, energy, process.cmME, weight, dim);

      // Set momenta for this event
      // for (int d = 0; d < dim; ++d) {
      //   for (int i = 0; i < 4; ++i) {
      //     for (int j = 0; j < 4; ++j) {
      //       printf(" p[%d][%d][%d] = %f",d,i,j,p(d,i,j));
      //     }
      //     printf("\n");
      //   }
      // }

      if (perf) {
        timer.Start();
      }
      // Evaluate matrix element
      // later process.sigmaKin(ncomb, goodhel, ntry, sum_hel, ngood, igood,
      // jhel);
      sigmaKin(p, meDevPtr, process.cHel, process.cIPD, process.cIPC, dim);//, debug, verbose);

      auto hp = Kokkos::create_mirror_view(p);
      Kokkos::deep_copy(hp,p);

      Kokkos::deep_copy(meHostPtr,meDevPtr);

      if (verbose)
        std::cout << "***********************************" << std::endl
                  << "Iteration #" << x+1 << " of " << numiter << std::endl;

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
                        << hp(d,i,0) << setiosflags(std::ios::scientific)
                        << std::setw(14) << hp(d,i,1)
                        << setiosflags(std::ios::scientific) << std::setw(14)
                        << hp(d,i,2) << setiosflags(std::ios::scientific)
                        << std::setw(14) << hp(d,i,3) << std::endl;
            std::cout << std::string(80, '-') << std::endl;
          }

          // Display matrix elements
          for (int i = 0; i < process.nprocesses; i++) {
            if (verbose)
              std::cout << " Matrix element = "
                        //	 << setiosflags(ios::fixed) << setprecision(17)
                        << meHostPtr(i*1 + d) << " GeV^" << meGeVexponent << std::endl;
            if (perf)
              matrixelementvector.push_back(meHostPtr(i*1 + d));
          }

          if (verbose)
            std::cout << std::string(80, '-') << std::endl;
        }
      } else if (!debug) {
        std::cout << ".";
      }

    } // end for numiter


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
      float sumelem = std::accumulate(matrixelementvector.begin(), matrixelementvector.end(), 0.0);
      float meanelem = sumelem / num_mes;
      float sqselem = std::inner_product(matrixelementvector.begin(), matrixelementvector.end(), 
                                         matrixelementvector.begin(), 0.0);
      float stdelem = std::sqrt(sqselem / num_mes - meanelem * meanelem);
      std::vector<double>::iterator maxelem = std::max_element(
          matrixelementvector.begin(), matrixelementvector.end());
      std::vector<double>::iterator minelem = std::min_element(
          matrixelementvector.begin(), matrixelementvector.end());

      std::cout << "***********************************" << std::endl
                << "NumIterations         = " << numiter << std::endl
                << "NumThreadsPerBlock    = " << gputhreads << std::endl
                << "NumBlocksPerGrid      = " << gpublocks << std::endl
                << "-----------------------------------" << std::endl
                << "NumberOfEntries       = " << num_wts << std::endl
                << std::scientific
                << "TotalTimeInWaveFuncs  = " << sum << " sec" << std::endl
                << "MeanTimeInWaveFuncs   = " << mean << " sec" << std::endl
                << "StdDevTimeInWaveFuncs = " << stdev << " sec" << std::endl
                << "MinTimeInWaveFuncs    = " << *mintime << " sec" << std::endl
                << "MaxTimeInWaveFuncs    = " << *maxtime << " sec" << std::endl
                << "-----------------------------------" << std::endl
                << "ProcessID:            = " << getpid() << std::endl
                << "NProcesses            = " << process.nprocesses << std::endl
                << "NumMatrixElements     = " << num_mes << std::endl
                << "MatrixElementsPerSec  = " << num_mes/sum << " sec^-1" << std::endl;

      std::cout << "***********************************" << std::endl
                << "NumMatrixElements     = " << num_mes << std::endl
                << std::scientific
                << "MeanMatrixElemValue   = " << meanelem << " GeV^" << meGeVexponent << std::endl
                << "StdErrMatrixElemValue = " << stdelem/sqrt(num_mes) << " GeV^" << meGeVexponent << std::endl
                << "StdDevMatrixElemValue = " << stdelem << " GeV^" << meGeVexponent << std::endl
                << "MinMatrixElemValue    = " << *minelem << " GeV^" << meGeVexponent << std::endl
                << "MaxMatrixElemValue    = " << *maxelem << " GeV^" << meGeVexponent << std::endl;
    }

  } // end Kokkos View Space
  Kokkos::finalize();
}
