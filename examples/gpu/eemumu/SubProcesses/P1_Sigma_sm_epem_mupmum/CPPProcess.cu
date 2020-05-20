#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include <algorithm>
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

__constant__ int cHel[16][4];
__constant__ double cMME[4];
__constant__ int cPerm[4];
__constant__ double cIPC[6];
__constant__ double cIPD[2];
// __shared__ double sw[6][12];

CPPProcess::CPPProcess(int numiterations, int gpublocks, int gputhreads,
                       bool verbose, bool debug)
    : m_numiterations(numiterations), gpu_nblocks(gpublocks),
      gpu_nthreads(gputhreads), m_verbose(verbose), m_debug(debug),
      dim(gpu_nblocks * gpu_nthreads), mME(4, 0.00) {

  static double tmpmME[4] = {0.0, 0.0, 0.0, 0.0};
  gpuErrchk(cudaMemcpyToSymbol(cMME, tmpmME, 4 * sizeof(double)));

  // Helicities for the process - nodim
  static int tHel[ncomb][nexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};
  gpuErrchk(cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int)));

  // perm - nodim
  static int perm[nexternal] = {0, 1, 2, 3};
  gpuErrchk(cudaMemcpyToSymbol(cPerm, perm, nexternal * sizeof(int)));
}

CPPProcess::~CPPProcess() {}

const std::vector<double> &CPPProcess::getMasses() const { return mME; }

void CPPProcess::initProc(std::string param_card_name) {

  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance();
  SLHAReader slha(param_card_name, m_verbose);
  pars->setIndependentParameters(slha);
  pars->setIndependentCouplings();
  static thrust::complex<double> tIPC[3] = {pars->GC_3, pars->GC_51,
                                            pars->GC_59};
  gpuErrchk(
      cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(thrust::complex<double>)));
  static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ};
  gpuErrchk(cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(double)));
  if (m_verbose) {
    pars->printIndependentParameters();
    pars->printIndependentCouplings();
  }
  // jamp2[0] = new double[1];
}

void CPPProcess::preSigmaKin() {
  // Set the parameters which change event by event
  pars->setDependentParameters();
  pars->setDependentCouplings();
  static bool firsttime = true;
  if (firsttime && m_verbose) {
    pars->printDependentParameters();
    pars->printDependentCouplings();
    firsttime = false;
  }
}
