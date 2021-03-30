#include <CL/sycl.hpp>
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
#include <cmath>

#define gpuErrchk3(ans)                                                        \
  { gpuAssert3((ans), __FILE__, __LINE__); }

class sycl_kernel;

inline void gpuAssert3(int code, const char *file, int line,
                       bool abort = true) {
}

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
  sycl::cpu_selector selector;
  sycl::queue q_ct1(selector);
  sycl::device dev_ct1 = q_ct1.get_device();;
  std::cout << "Device: " << dev_ct1.get_info<sycl::info::device::name>() << "\n";
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

  if (verbose)
    std::cout << "# iterations: " << numiter << std::endl;

  // Create a process object
  CPPProcess process(numiter, gpublocks, gputhreads, verbose, debug);

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  double energy = 1500;
  double weight;

  int meGeVexponent = -(2 * process.nexternal - 8);

  int dim = gpublocks * gputhreads;

  // Local Memory
  //typedef double arr_t[6][4];
  double* lp = new double[6*3*dim];

  std::vector<double> matrixelementvector;

	  
  
  for (int x = 0; x < numiter; ++x) {
    // Get phase space point
    std::vector<std::vector<double *> > p = get_momenta(process.ninitial, energy, process.getMasses(), weight, dim);

    // Set momenta for this event
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
          lp[i*dim*3+j*dim+d] = p[d][i][1+j];
        }
      }
    }
    //new
    double * allmomenta = sycl::malloc_shared<double>(3*6*dim*sizeof(double) ,q_ct1);
    double * me = sycl::malloc_shared<double>(dim*sizeof(double), q_ct1);

  if (perf) {
	  timer.Start();
  }
  
  // Evaluate matrix element
  // later process.sigmaKin(ncomb, goodhel, ntry, sum_hel, ngood, igood,
  // jhel);
  
  /* Setup Kernel Variables */
  sycl::nd_range nd_range = sycl::nd_range<3>(sycl::range<3>(1, 1, gpublocks) *
                                              sycl::range<3>(1, 1, gputhreads),
                                              sycl::range<3>(1, 1, gputhreads));
  
  // USM Allocations, copy host -> device
  const int cHEL_ROW = 64;
  const int cHEL_COL = 6;
  int **cHel_local = sycl::malloc_device<int*>(cHEL_ROW *sizeof(int*), q_ct1);
  // allocate elements using host pointers, copy
  for (int i=0; i< cHEL_ROW ; i++){
	  cHel_local[i] = sycl::malloc_device<int>(cHEL_COL * sizeof(int), q_ct1);
	  q_ct1.memcpy(cHel_local[i], &cHel[i][0], cHEL_COL * sizeof(int)).wait_and_throw();	  
  }
  // copy host-side pointers
  q_ct1.memcpy(cHel_local, &cHel[0], cHEL_ROW * cHEL_COL * sizeof(int)).wait_and_throw();

  double *cIPC_local = sycl::malloc_device<double>(sizeof(cIPC), q_ct1);
  q_ct1.memcpy(cIPC_local, cIPC, sizeof(cIPC)).wait_and_throw();
        
  double *cIPD_local = sycl::malloc_device<double>(sizeof(cIPD), q_ct1);
  q_ct1.memcpy(cIPD_local, cIPD, sizeof(cIPD)).wait_and_throw();

  try{
	  q_ct1.submit([&](sycl::handler &cgh) {
			  // instantiate kernel
			  auto k = sigmaKin(allmomenta,
			                    me,
			                    cHel_local,
			                    cIPC_local,
			                    cIPD_local
			                    );
			  
			  cgh.parallel_for(nd_range, k);

		  }).wait_and_throw(); //, debug, verbose);
  } catch(const sycl::exception& e){
	  std::cout << "error: " << e.what();
  }

  // Shared matrix element means no need to memcpy to host

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
                      << me[i*1 + d] << " GeV^" << meGeVexponent << std::endl;
          if (perf)
            matrixelementvector.push_back(me[i*1 + d]);
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }
    } else if (!debug) {
      std::cout << ".";
    }

    for (std::vector<std::vector<double *> >::iterator it = p.begin();
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
              << std::scientific << "MeanMatrixElemValue   = " << meanelem
              << " GeV^" << meGeVexponent << std::endl
              << "StdErrMatrixElemValue = " << stdelem / sqrt(num_mes)
              << " GeV^" << meGeVexponent << std::endl
              << "StdDevMatrixElemValue = " << stdelem << " GeV^"
              << meGeVexponent << std::endl
              << "MinMatrixElemValue    = " << *minelem << " GeV^"
              << meGeVexponent << std::endl
              << "MaxMatrixElemValue    = " << *maxelem << " GeV^"
              << meGeVexponent << std::endl;
  }
  delete[] lp;

}
