//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.0, 2020-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
// #include "HelAmps_sm.h"
#include "HelAmps_sm.h"

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

  // mm
  // mME
  static double tmpmME[4] = {0.00, 0.00, 0.00, 0.00};
  gpuErrchk(cudaMallocManaged(&m->mME, 4 * sizeof(double)));
  gpuErrchk(cudaMemcpy((void *)m->mME, (void *)tmpmME, 4 * sizeof(double),
                       cudaMemcpyHostToHost));

  // mm
  // p
  ptemp = (double **)malloc(m->tnioparticles * sizeof(double));
  for (int i = 0; i < 4; ++i) {
    ptemp[i] = (double *)malloc(4 * sizeof(double));
  }
  gpuErrchk(cudaMallocManaged(&m->tp, m->tnioparticles * sizeof(double)));
  for (int i = 0; i < 4; ++i) {
    gpuErrchk(cudaMallocManaged(&m->tp[i], 4 * sizeof(double)));
  }

  // amp
  gpuErrchk(cudaMallocManaged(&m->tamp, m->tnamplitudes *
                                            sizeof(thrust::complex<double>)));

  // w
  gpuErrchk(
      cudaMallocManaged(&m->tw, m->twrows * sizeof(thrust::complex<double>)));
  for (int i = 0; i < m->twrows; ++i) {
    gpuErrchk(cudaMallocManaged(
        &m->tw[i], m->tnwavefuncs * sizeof(thrust::complex<double>)));
    for (int y = 0; y < m->twrows; ++y) {
      m->tw[i][y] = thrust::complex<double>(0.00, 0.00);
    }
  }

  // mm
  // helicities
  // Helicities for the process
  static int helicities[m->tncomb][m->tnexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};
  gpuErrchk(cudaMallocManaged(&m->thelicities,
                              m->tncomb * m->tnexternal * sizeof(int)));
  for (int i = 0; i < m->tncomb; ++i) {
    gpuErrchk(cudaMemcpy(&m->thelicities[i], &helicities[i],
                         m->tnexternal * sizeof(int), cudaMemcpyHostToHost));
  }

  // mm
  // perm
  static int perm[m->tnexternal];
  for (int i = 0; i < m->tnexternal; i++) {
    perm[i] = i;
  }
  gpuErrchk(cudaMallocManaged(&m->tperm, m->tnexternal * sizeof(int)));
  gpuErrchk(cudaMemcpy((void *)m->tperm, (void *)perm,
                       m->tnexternal * sizeof(int), cudaMemcpyHostToHost));
}

CPPProcess::~CPPProcess() {}

void CPPProcess::setMomenta(std::vector<double *> &momenta) {

  for (int i = 0; i < m->tnioparticles; ++i) {
    gpuErrchk(cudaMemcpy((void *)m->tp[i], (void *)momenta[i],
                         4 * sizeof(double), cudaMemcpyHostToHost));
    // memcpy(ptemp[i], momenta[i], 4 * sizeof(double));
  }

  /*
  for (std::vector<double *>::iterator it = momenta.begin();
       it != momenta.end(); ++it) {
    double *tmp;
    cudaMallocManaged(&tmp, 4 * sizeof(double));
    cudaMemcpy(tmp, *it, 4 * sizeof(double), cudaMemcpyHostToHost);
    p.push_back(tmp);
  }
  */
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

void CPPProcess::resetGPUMemory() {

  for (std::vector<double *>::iterator it = p.begin(); it != p.end(); ++it) {
    gpuErrchk(cudaFree(*it));
  }
  p.clear();
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

  gMG5_sm::calculate_wavefunctions<<<1, 1>>>(
      m->tperm, m->thelicities[ihel], m->mME, m->tp, m->tamp, m->tw, pars->GC_3,
      pars->GC_51, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ);
  cudaDeviceSynchronize();

  gpuErrchk(cudaMemcpy((void *)amp, (void *)m->tamp,
                       m->tnamplitudes * sizeof(thrust::complex<double>),
                       cudaMemcpyDeviceToHost));

  float gputime = m_timer.GetDuration();
  std::cout << "Wave function time: " << gputime << std::endl;
}

/*
void CPPProcess::calculate_wavefunctions(const int perm[], const int hel[]) {
  // Calculate wavefunctions for all processes
  // int i, j;


  // Calculate all wavefunctions
  gMG5_sm::oxxxxx<<<1, 1>>>(p[perm[0]], mME[0], hel[0], -1, 0);
  gMG5_sm::ixxxxx<<<1, 1>>>(p[perm[1]], mME[1], hel[1], +1, 1);
  gMG5_sm::ixxxxx<<<1, 1>>>(p[perm[2]], mME[2], hel[2], -1, 2);
  gMG5_sm::oxxxxx<<<1, 1>>>(p[perm[3]], mME[3], hel[3], +1, 3);

  gMG5_sm::FFV1P0_3<<<1, 1>>>(1, 0, pars->GC_3, pars->ZERO, pars->ZERO, 4);
  gMG5_sm::FFV2_4_3<<<1, 1>>>(1, 0, -pars->GC_51, pars->GC_59, pars->mdl_MZ,
                              pars->mdl_WZ, 5);
  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  gMG5_sm::FFV1_0<<<1, 1>>>(2, 3, 4, pars->GC_3, &m->tamp[0]);
  gMG5_sm::FFV2_4_0<<<1, 1>>>(2, 3, 5, -pars->GC_51, pars->GC_59, &m->tamp[1]);
  cudaDeviceSynchronize();

}
*/

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
  jamp[0] = -m->tamp[0] - m->tamp[1];

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
