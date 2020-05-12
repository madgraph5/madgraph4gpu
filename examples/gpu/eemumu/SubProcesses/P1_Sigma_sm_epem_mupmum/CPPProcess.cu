//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include <algorithm>
#include <algorithm> // perf stats
#include <iomanip>   // setw
#include <iostream>
#include <numeric> // perf stats
#include <thrust/complex.h>

#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

__constant__ int hel2[16][4];

CPPProcess::CPPProcess(int numiterations, int gpuwarps, int gputhreads,
                       bool verbose, bool debug, bool perf)
    : m_numiterations(numiterations), gpu_nwarps(gpuwarps),
      gpu_nthreads(gputhreads), m_verbose(verbose), m_debug(debug),
      m_perf(perf), dim(gpu_nwarps * gpu_nthreads), mME(4, 0.00) {

  amp = new thrust::complex<double> *[dim];
  for (int i = 0; i < dim; ++i) {
    amp[i] = new thrust::complex<double>[namplitudes];
  }

  matrix_element = new double *[dim];
  for (int i = 0; i < dim; ++i) {
    matrix_element[i] = new double[nprocesses];
  }

  m = new processMem();

  static double tmpmME[4] = {0.00, 0.00, 0.00, 0.00};
  gpuErrchk(cudaMalloc(&m->tmME, 4 * sizeof(double)));
  gpuErrchk(cudaMemcpy((void *)m->tmME, (void *)tmpmME, 4 * sizeof(double),
                       cudaMemcpyHostToDevice));

  /*
  - rambo::get_momenta (rambo::rambo) fills vector 4 particles * 4 momenta
  - CPPProcess::setMomenta fills m->tp
  */
  gpuErrchk(cudaMalloc(&m->tp, dim * nioparticles * 4 * sizeof(double)));

  /*
  - i/o parameter to kernels
  - after kernel call fill amp member
  - matrix_1_epem_mupmum uses amp for final calculation of weight -> return
  value used in CPPProcess::sigmaKin --> writes into CPPProcess::matrix_element
  member
  */
  gpuErrchk(cudaMallocManaged(&m->tamp, dim * namplitudes *
                                            sizeof(thrust::complex<double>)));

  /*
  temporary variables needed inside the kernel wave functions (re-using them
  inside two functions )
  */
  /*
  gpuErrchk(cudaMalloc(&m->tmp, dim * 4 * sizeof(thrust::complex<double>)));
  */

  // Helicities for the process - nodim
  static int helicities[ncomb][nexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};
  gpuErrchk(cudaMemcpyToSymbol(hel2, helicities, 16 * 4 * sizeof(int)));

  // perm - nodim
  static int perm[nexternal];
  for (int i = 0; i < nexternal; i++)
    perm[i] = i;
  gpuErrchk(cudaMalloc(&m->tperm, nexternal * sizeof(int)));
  gpuErrchk(cudaMemcpy((void *)m->tperm, (void *)perm, nexternal * sizeof(int),
                       cudaMemcpyHostToDevice));
}

CPPProcess::~CPPProcess() {}

void CPPProcess::setMomenta(std::vector<std::vector<double *>> &momenta) {

  for (int d = 0; d < dim; ++d) {
    for (int i = 0; i < nioparticles; ++i) {
      // cudaDeviceSynchronize();
      gpuErrchk(cudaMemcpy((void *)m->tp[d][i], (void *)momenta[d][i],
                           4 * sizeof(double), cudaMemcpyHostToDevice));
      // cudaDeviceSynchronize();
    }
  }
}

const std::vector<double> &CPPProcess::getMasses() const { return mME; }

void CPPProcess::printPerformanceStats() {
  float sum = std::accumulate(m_wavetimes.begin(), m_wavetimes.end(), 0.0);
  int numelems = m_wavetimes.size();
  float mean = sum / numelems;
  float sq_sum = std::inner_product(m_wavetimes.begin(), m_wavetimes.end(),
                                    m_wavetimes.begin(), 0.0);
  float stdev = std::sqrt(sq_sum / numelems - mean * mean);
  std::vector<float>::iterator mintime =
      std::min_element(m_wavetimes.begin(), m_wavetimes.end());
  std::vector<float>::iterator maxtime =
      std::max_element(m_wavetimes.begin(), m_wavetimes.end());

  std::cout << "***********************************" << std::endl
            << "NumIterations        = " << m_numiterations << std::endl
            << "NumThreadsPerBlock   = " << gpu_nthreads << std::endl
            << "NumBlocksPerGrid     = " << gpu_nwarps << std::endl
            << "NumberOfEntries      = " << numelems << std::endl
            << std::scientific << "TotalTimeInWaveFuncs = " << sum << std::endl
            << "MeanTimeinWaveFuncs  = " << mean << std::endl
            << "StdDevWaveFuncs      = " << stdev << std::endl
            << "MinTimeInWaveFuncs   = " << *mintime << std::endl
            << "MaxTimeInWaveFuncs   = " << *maxtime << std::endl;
}

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(std::string param_card_name) {

  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance();
  SLHAReader slha(param_card_name, m_verbose);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
  if (m_verbose) {
    pars->printIndependentParameters();
    pars->printIndependentCouplings();
  }
  jamp2[0] = new double[1];
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess::sigmaKin() {
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
  static bool firsttime = true;
  if (firsttime && m_verbose) {
    pars->printDependentParameters();
    pars->printDependentCouplings();
    firsttime = false;
  }

  // Reset color flows
  for (int i = 0; i < 1; i++)
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 16;
  static bool goodhel[ncomb] = {ncomb * false};
  static int ntry = 0, sum_hel = 0, ngood = 0;
  static int igood[ncomb];
  static int jhel;
  // sr fixme // move out of functions
  double **t = new double *[nprocesses];
  for (int i = 0; i < nprocesses; ++i) {
    t[i] = new double[dim];
  }
  // double t[nprocesses][dim] = {0}; // [nprocesses]
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {4};

  ntry = ntry + 1;

  // Reset the matrix elements
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < nprocesses; ++j) {
      matrix_element[i][j] = 0.;
    }
  }

  // sr fixme // better to run the first n calculations serial?
  if (sum_hel == 0 || ntry < 10) {
    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
      if (goodhel[ihel] || ntry < 2) {

        call_wavefunctions_kernel(ihel);
        matrix_1_epem_mupmum(t[0]);

        double tsum[dim] = {0};
        for (int iproc = 0; iproc < nprocesses; iproc++) {
          for (int d = 0; d < dim; ++d) {
            matrix_element[d][iproc] += t[iproc][d];
            tsum[d] += t[iproc][d];
          }
        }
        // Store which helicities give non-zero result
        for (int d = 0; d < dim; ++d) {
          if (tsum[d] != 0. && !goodhel[ihel]) {
            goodhel[ihel] = true;
            ngood++;
            igood[ngood] = ihel;
          }
        }
      }
    }
    jhel = 0;
    sum_hel = min(sum_hel, ngood);
  } else {
    // Only use the "good" helicities
    // sr fixme // is the calculation of good helicities parralelizable?
    for (int j = 0; j < sum_hel; j++) {
      jhel++;
      if (jhel >= ngood)
        jhel = 0;
      double hwgt = double(ngood) / double(sum_hel);
      int ihel = igood[jhel];

      call_wavefunctions_kernel(ihel);
      matrix_1_epem_mupmum(t[0]);

      for (int iproc = 0; iproc < nprocesses; iproc++) {
        for (int d = 0; d < dim; ++d) {
          matrix_element[d][iproc] += t[iproc][d] * hwgt;
        }
      }
    }
    // sr fixme // move out of functions
    for (int i = 0; i < nprocesses; ++i) {
      delete[] t[i];
    }
    delete[](*t);
  }

  for (int i = 0; i < dim; i++)
    for (int j = 0; j < nprocesses; ++j) {
      matrix_element[i][j] /= denominators[j];
    }
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess::call_wavefunctions_kernel(int ihel) {

  if (m_perf) {
    m_timer.Start();
  }

  // cudaDeviceSynchronize();
  gMG5_sm::calculate_wavefunctions<<<gpu_nwarps, gpu_nthreads>>>(
      m->tperm, ihel, m->tmME, m->tp, m->tamp, pars->GC_3, pars->GC_51,
      pars->GC_59, pars->mdl_MZ, pars->mdl_WZ, m_debug, m_verbose);
  cudaDeviceSynchronize();

  // memcpy(amp, m->tamp, namplitudes * sizeof(thrust::complex<double>));
  /*
  gpuErrchk(cudaMemcpy(amp, m->tamp,
                       dim * namplitudes * sizeof(thrust::complex<double>),
                       cudaMemcpyDeviceToHost));
*/
  if (m_perf) {
    float gputime = m_timer.GetDuration();
    m_wavetimes.push_back(gputime);
    if (m_verbose)
      std::cout << "Wave function time: " << gputime << std::endl;
  }

  /*
    std::cout << std::endl << std::endl;
    for (int i = 0; i < gpu_nwarps; ++i) {
      for (int j = 0; j < gpu_nthreads; ++j) {
        int d = i * gpu_nwarps + j;
        std::cout << m->tamp[d][0] << ", " << m->tamp[d][1] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  */
}

// --> calculate multi-dimensional amp
void CPPProcess::matrix_1_epem_mupmum(double *matrix) {
  // static double matrix[dim];
  for (int d = 0; d < dim; ++d) {
    thrust::complex<double> *lamp = m->tamp[d]; // [namplitudes]
    int i, j;
    // Local variables
    // const int ngraphs = 2;
    const int ncolor = 1;
    std::complex<double> ztemp;
    std::complex<double> jamp[ncolor];
    // The color matrix;
    static const double denom[ncolor] = {1};
    static const double cf[ncolor][ncolor] = {{1}};

    // Calculate color flows
    jamp[0] = -lamp[0] - lamp[1];

    // Sum and square the color flows to get the matrix element
    matrix[d] = 0;
    for (i = 0; i < ncolor; i++) {
      ztemp = 0.;
      for (j = 0; j < ncolor; j++)
        ztemp = ztemp + cf[i][j] * jamp[j];
      matrix[d] = matrix[d] + real(ztemp * conj(jamp[i])) / denom[i];
    }

    // Store the leading color flows for choice of color
    // sr fixme // maybe this needs to go outside the loop? does it need a
    // dimension?
    for (i = 0; i < ncolor; i++)
      jamp2[0][i] += real(jamp[i] * conj(jamp[i]));
  }
}
