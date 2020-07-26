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

int usage(char* argv0, int ret = 1) {
  std::cout << "Usage: " << argv0 
            << " [--verbose|-v] [--debug|-d] [--performance|-p]"
            << " [#gpuBlocksPerGrid #gpuThreadsPerBlock] #iterations" << std::endl;
  return ret;
}

int main(int argc, char **argv) {
  bool verbose = false, debug = false, perf = false;
  int niter = 0;
  int gpublocks = 1;
  int gputhreads = 1;
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
    niter = numvec[2];
  } else if (veclen == 1) {
    niter = numvec[0];
  } else {
    return usage(argv[0]);
  }

  if (niter == 0)
    return usage(argv[0]);

  //std::cout << "Calling cudaFree... " << std::endl;
  cudaFree(0); // SLOW!
  //std::cout << "Calling cudaFree... done" << std::endl;

  if (verbose)
    std::cout << "# iterations: " << niter << std::endl;

  // Create a process object
  CPPProcess process( niter, gpublocks, gputhreads, verbose, debug );

  // Read param_card and set parameters
  process.initProc("../../Cards/param_card.dat");

  double energy = 1500;
  double weight;

  int meGeVexponent = -(2 * process.nexternal - 8);

  int ndim = gpublocks * gputhreads;

  // Local Memory
  const int npar = process.nexternal; // for this process (eemumu): npar=4
  const int np3 = 3; // dimension of 3-momenta (px,py,pz): energy from rambo is ignored (mass is known)
  double* hstMomenta = new double[npar*np3*ndim];

  double* meHostPtr = new double[ndim*1];
  double *meDevPtr =0;
  int num_bytes_back = 1 * ndim * sizeof(double);
  cudaMalloc((void**)&meDevPtr, num_bytes_back);


  std::vector<double> matrixelementvector;


  for ( int iiter = 0; iiter < niter; ++iiter ) 
  {

    //std::cout << "Iteration #" << iiter+1 << " of " << niter << std::endl;
    // Get a vector of ndim phase space points
    std::vector<std::vector<double *>> rmbMomenta = // AOS[ndim][npar][np3+1]
      get_momenta(process.ninitial, energy, process.getMasses(), weight, ndim); // SLOW!
    //std::cout << "Got momenta" << std::endl;

    // Set momenta for this event by copying them from the rambo output
    for (int idim = 0; idim < ndim; ++idim)
      for (int ipar = 0; ipar < npar; ++ipar)
        for (int ip3 = 0; ip3 < np3; ++ip3)
          hstMomenta[ipar*ndim*np3 + ip3*ndim + idim] = // SOA[npar][np3][ndim]
            rmbMomenta[idim][ipar][1+ip3]; // energy from rambo is ignored (mass is known)

    int nbytesMomenta = np3*npar*ndim * sizeof(double);
    double* devMomenta = 0; // device momenta
    cudaMalloc( &devMomenta, nbytesMomenta );
    cudaMemcpy( devMomenta, hstMomenta, nbytesMomenta, cudaMemcpyHostToDevice );

    //gpuErrchk3(cudaMemcpy3D(&tdp));

   //process.preSigmaKin();

    if (perf) {
      timer.Start();
    }

    // Evaluate matrix element
    // later process.sigmaKin(ncomb, goodhel, ntry, sum_hel, ngood, igood,
    // jhel);
    sigmaKin<<<gpublocks, gputhreads>>>(devMomenta,  meDevPtr);//, debug, verbose);
    gpuErrchk3( cudaPeekAtLastError() );
    //gpuErrchk3(cudaMemcpy2D(meHostPtr, sizeof(double), meDevPtr, mePitch,
    //                        sizeof(double), dim, cudaMemcpyDeviceToHost));

   cudaMemcpy(meHostPtr, meDevPtr, 1 * ndim*sizeof(double), cudaMemcpyDeviceToHost);

    if (verbose)
      std::cout << "***********************************" << std::endl
                << "Iteration #" << iiter+1 << " of " << niter << std::endl;

    if (perf) {
      float gputime = timer.GetDuration();
      wavetimes.push_back(gputime);
      if (verbose)
        std::cout << "Wave function time: " << gputime << std::endl;
    }

    if (verbose || perf) {

      for (int idim = 0; idim < ndim; ++idim) {

        if (verbose) {
          std::cout << "Momenta:" << std::endl;
          for (int ipar = 0; ipar < npar; ipar++)
            std::cout << std::setw(4) << ipar + 1
                      << setiosflags(std::ios::scientific)
                      << std::setw(14) << rmbMomenta[idim][ipar][0]
                      << setiosflags(std::ios::scientific)
                      << std::setw(14) << rmbMomenta[idim][ipar][1]
                      << setiosflags(std::ios::scientific)
                      << std::setw(14) << rmbMomenta[idim][ipar][2]
                      << setiosflags(std::ios::scientific)
                      << std::setw(14) << rmbMomenta[idim][ipar][3] << std::endl;
          std::cout << std::string(80, '-') << std::endl;
        }

        // Display matrix elements
        for (int iproc = 0; iproc < process.nprocesses; iproc++) {
          if (verbose)
            std::cout << " Matrix element = "
                      //	 << setiosflags(ios::fixed) << setprecision(17)
                      << meHostPtr[iproc*1 + idim] << " GeV^" << meGeVexponent << std::endl;
          if (perf)
            matrixelementvector.push_back(meHostPtr[iproc*1 + idim]);
        }

        if (verbose)
          std::cout << std::string(80, '-') << std::endl;
      }
    } else if (!debug) {
      std::cout << ".";
    }

    for (std::vector<std::vector<double *>>::iterator it = rmbMomenta.begin();
         it != rmbMomenta.end(); ++it) {
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
              << "NumIterations         = " << niter << std::endl
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
  delete[] hstMomenta;

}
