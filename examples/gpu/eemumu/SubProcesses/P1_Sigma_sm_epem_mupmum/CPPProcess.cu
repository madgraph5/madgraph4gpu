//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include <algorithm>
#include <iomanip> // setw
#include <iostream>
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

CPPProcess::CPPProcess(bool verbose, bool debug)
    : m_verbose(verbose), m_debug(debug), mME(4, 0.00) {

  m = new processMem();
  int dim = gpu_nwarps * gpu_nthreads;

  // tmME - nodim
  static double tmpmME[4] = {0.00, 0.00, 0.00, 0.00};
  gpuErrchk(cudaMalloc(&m->tmME, 4 * sizeof(double)));
  gpuErrchk(cudaMemcpy((void *)m->tmME, (void *)tmpmME, 4 * sizeof(double),
                       cudaMemcpyHostToDevice));

  // tp - DIM --> LATER --> add "dim" to cudaMalloc
  /*
  - rambo::get_momenta (rambo::rambo) fuellt vector 4 particles * 4 momenta
  - CPPProcess::setMomenta fuellt m->tp
  */
  gpuErrchk(cudaMalloc(&m->tp, nioparticles * 4 * sizeof(double)));

  // amp - DIM --> LATER --> add "dim" to cudaMalloc
  /*
  - i/o parameter to kernels
  - after kernel call fill amp member
  - matrix_1_epem_mupmum uses amp for final calculation of weight -> return
  value used in rambo::sigmaKin --> writes into CPPProcess::matrix_element
  member
  */
  gpuErrchk(
      cudaMalloc(&m->tamp, namplitudes * sizeof(thrust::complex<double>)));

  // w - DIM
  /*
  - internal variable within kernels
  */
  gpuErrchk(cudaMalloc((void **)&m->tw, dim * wrows * nwavefuncs *
                                            sizeof(thrust::complex<double>)));

  // Helicities for the process - nodim
  static int helicities[ncomb][nexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};
  gpuErrchk(cudaMalloc(&m->thelicities, ncomb * nexternal * sizeof(int)));
  gpuErrchk(cudaMemcpy(m->thelicities, helicities,
                       ncomb * nexternal * sizeof(int),
                       cudaMemcpyHostToDevice));

  // perm - nodim
  static int perm[nexternal];
  for (int i = 0; i < nexternal; i++)
    perm[i] = i;
  gpuErrchk(cudaMalloc(&m->tperm, nexternal * sizeof(int)));
  gpuErrchk(cudaMemcpy((void *)m->tperm, (void *)perm, nexternal * sizeof(int),
                       cudaMemcpyHostToDevice));
}

CPPProcess::~CPPProcess() {}

void CPPProcess::setMomenta(std::vector<double *> &momenta) {

  for (int i = 0; i < nioparticles; ++i) {
    gpuErrchk(cudaMemcpy((void *)m->tp[i], (void *)momenta[i],
                         4 * sizeof(double), cudaMemcpyHostToDevice));
  }
}

const std::vector<double> &CPPProcess::getMasses() const { return mME; }

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
  // thrust::complex<double> **wfs;
  double t[nprocesses];
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {4};

  ntry = ntry + 1;

  // Reset the matrix elements
  for (int i = 0; i < nprocesses; i++) {
    matrix_element[i] = 0.;
  }

  if (sum_hel == 0 || ntry < 10) {
    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
      if (goodhel[ihel] || ntry < 2) {

        call_wavefunctions_kernel(ihel);
        t[0] = matrix_1_epem_mupmum();

        double tsum = 0;
        for (int iproc = 0; iproc < nprocesses; iproc++) {
          matrix_element[iproc] += t[iproc];
          tsum += t[iproc];
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel]) {
          goodhel[ihel] = true;
          ngood++;
          igood[ngood] = ihel;
        }
      }
    }
    jhel = 0;
    sum_hel = min(sum_hel, ngood);
  } else {
    // Only use the "good" helicities
    for (int j = 0; j < sum_hel; j++) {
      jhel++;
      if (jhel >= ngood)
        jhel = 0;
      double hwgt = double(ngood) / double(sum_hel);
      int ihel = igood[jhel];

      call_wavefunctions_kernel(ihel);
      t[0] = matrix_1_epem_mupmum();

      for (int iproc = 0; iproc < nprocesses; iproc++) {
        matrix_element[iproc] += t[iproc] * hwgt;
      }
    }
  }

  for (int i = 0; i < nprocesses; i++)
    matrix_element[i] /= denominators[i];
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess::sigmaHat() {
  // Select between the different processes
  if (id1 == -11 && id2 == 11) {
    // Add matrix elements for processes with beams (-11, 11)
    return matrix_element[0];
  } else {
    // Return 0 if not correct initial state assignment
    return 0.;
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess::call_wavefunctions_kernel(int ihel) {

  m_timer.Start();

  gMG5_sm::calculate_wavefunctions<<<gpu_nwarps, gpu_nthreads>>>(
      m->tperm, m->thelicities, ihel, m->tmME, m->tp, m->tamp, m->tw,
      pars->GC_3, pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ);
  // cudaDeviceSynchronize();

  // memcpy(amp, m->tamp, namplitudes * sizeof(thrust::complex<double>));
  gpuErrchk(cudaMemcpy((void *)amp, (void *)m->tamp,
                       namplitudes * sizeof(thrust::complex<double>),
                       cudaMemcpyDeviceToHost));

  float gputime = m_timer.GetDuration();
  std::cout << "Wave function time: " << gputime << std::endl;
}

double CPPProcess::matrix_1_epem_mupmum() {
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
  jamp[0] = -amp[0] - amp[1];

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for (i = 0; i < ncolor; i++) {
    ztemp = 0.;
    for (j = 0; j < ncolor; j++)
      ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + real(ztemp * conj(jamp[i])) / denom[i];
  }

  // Store the leading color flows for choice of color
  for (i = 0; i < ncolor; i++)
    jamp2[0][i] += real(jamp[i] * conj(jamp[i]));

  return matrix;
}
