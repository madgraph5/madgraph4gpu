#include <algorithm> // perf stats
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric> // perf stats
#include <unistd.h>
#include <vector>

#include <CL/sycl.hpp>
#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include "rambo.h"
#include "timer.h"

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
            << " [#work items] #iterations" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false, debug = false, perf = false;
  int numiter = 0, work_groups = 0, work_items = 1;
  std::vector<int> numvec;
  Timer<TIMERTYPE> timer;
  std::vector<float> wavetimes;

  // remove
  // int nprocesses = 1;

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
  if (veclen == 2) {
    work_items = numvec[0];
    numiter = numvec[1];
  } else if (veclen == 1) {
    numiter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (numiter == 0)
    return usage(argv[0]);

  cl::sycl::queue q; 
  cl::sycl::device device = q.get_device();
  int max_cus = device.get_info<cl::sycl::info::device::max_compute_units>();
  if (work_items >  max_cus){
      std::cerr << "Error: Specified work items is greater than available: " << work_items 
        << " < " << max_cus << "."<< std::endl;
        return 1;
  }

  if (verbose)
    std::cout << "# iterations: " << numiter << std::endl;

  // Create a process object
  CPPProcess process(numiter, work_groups, work_items, verbose, debug);

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  double energy = 1500;
  double weight;

  int meGeVexponent = -(2 * process.nexternal - 8);

  int dim = work_items;

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

    if (perf) {
      timer.Start();
    }

    cl::sycl::range<1> range{(unsigned long)dim};

    double* matrix_element[dim];
 
    q.submit([&](cl::sycl::handler& cgh){

      // Put here so that it captures the instance
      double lp[1][4][4];
      // Set momenta for this event
      for (int d = 0; d < 1; ++d) {
        for (int i = 0; i < 4; ++i) {
          for (int j = 0; j < 4; ++j) {
            lp[d][i][j] = p[d][i][j];
          }
        }
      }
      std::complex<double> IPC[3];
      IPC[0] = process.IPC[0];
      IPC[1] = process.IPC[1];
      IPC[2] = process.IPC[2];
      double IPD[2];
      IPD[0] = process.IPD[0];
      IPD[1] = process.IPD[1]; 

      process.preSigmaKin();
    cgh.parallel_for<class my_kernel>(range,
                                         [=] (cl::sycl::id<1> idx) {
       gMG5_sm::sigmaKin(lp[0], IPC, IPD, debug, verbose);
        
    
      }); // End parallel_for
    }); // End submit
    q.wait();
    
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
                      << matrix_element[d][i] << " GeV^" << meGeVexponent << std::endl;
          //if (perf)
          //  matrixelementvector.push_back(meHostPtr[d][i]);
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }
    } else if (!debug) {
     // std::cout << ".";
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
              << "NumThreadsPerBlock    = " << work_items << std::endl
              << "NumBlocksPerGrid      = " << work_groups << std::endl
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
}
