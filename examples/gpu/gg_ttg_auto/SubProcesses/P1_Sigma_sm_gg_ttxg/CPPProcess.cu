//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "../../src/HelAmps_sm.h"

#include <algorithm> 
#include <iostream> 
#include <thrust/complex.h> 

#include "../../src/HelAmps_sm.cu"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ g WEIGHTED<=3 @1

__constant__ int cHel[32][5]; 
// __constant__ double cmME[5]; value hardcoded now
// extern __constant__ int cPerm[4];
// 
__constant__ double cIPC[6];  // coupling ?
__constant__ double cIPD[2]; 


// Evaluate |M|^2 for each subprocess

__device__ void calculate_wavefunctions(int ihel, double local_mom[5][3],
    double &matrix)
{
  thrust::complex<double> amp[1];  // was 18
  const int ncolor = 6; 
  thrust::complex<double> jamp[ncolor]; 
  // Calculate wavefunctions for all processes
  using namespace MG5_sm; 
  thrust::complex<double> w[12][6]; 
  for(int i = 0; i < 6; i++ )
  {
    jamp[i] = thrust::complex<double> (0., 0.); 
  }
  vxxxxx(local_mom[0], 0., cHel[ihel][0], -1, w[0]); 
  vxxxxx(local_mom[1], 0., cHel[ihel][1], -1, w[1]); 
  oxxxxx(local_mom[2], cIPD[0], cHel[ihel][2], +1, w[2]); 
  ixxxxx(local_mom[3], cIPD[0], cHel[ihel][3], -1, w[3]); 
  vxxxxx(local_mom[4], 0., cHel[ihel][4], +1, w[4]); 
  VVV1P0_1(w[0], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[5]);
  FFV1P0_3(w[3], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[6]);
  // Amplitude(s) for diagram number 1
  VVV1_0(w[5], w[6], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[0]);
  jamp[0] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[4] += +amp[0]; 
  jamp[5] += -amp[0]; 
  FFV1_1(w[2], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[7]);
  // Amplitude(s) for diagram number 2
  FFV1_0(w[3], w[7], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[4] += +thrust::complex<double> (0, 1) * amp[0]; 
  jamp[5] += -thrust::complex<double> (0, 1) * amp[0]; 
  FFV1_2(w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[8]);
  // Amplitude(s) for diagram number 3
  FFV1_0(w[8], w[2], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[0] += +thrust::complex<double> (0, 1) * amp[0]; 
  jamp[2] += -thrust::complex<double> (0, 1) * amp[0]; 
  FFV1_1(w[2], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[5]);
  FFV1_2(w[3], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[9]);
  // Amplitude(s) for diagram number 4
  FFV1_0(w[9], w[5], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[1] += -amp[0]; 
  VVV1P0_1(w[1], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[10]);
  // Amplitude(s) for diagram number 5
  FFV1_0(w[3], w[5], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[0] += +thrust::complex<double> (0, 1) * amp[0]; 
  jamp[1] += -thrust::complex<double> (0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 6
  FFV1_0(w[8], w[5], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[0] += -amp[0]; 
  FFV1_2(w[3], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[5]);
  FFV1_1(w[2], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[11]);
  // Amplitude(s) for diagram number 7
  FFV1_0(w[5], w[11], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[3] += -amp[0]; 
  // Amplitude(s) for diagram number 8
  FFV1_0(w[5], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[3] += +thrust::complex<double> (0, 1) * amp[0]; 
  jamp[5] += -thrust::complex<double> (0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 9
  FFV1_0(w[5], w[7], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[5] += -amp[0]; 
  VVV1P0_1(w[0], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[5]);
  // Amplitude(s) for diagram number 10
  FFV1_0(w[3], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[2] += +thrust::complex<double> (0, 1) * amp[0]; 
  jamp[3] += -thrust::complex<double> (0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 11
  FFV1_0(w[9], w[2], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[1] += +thrust::complex<double> (0, 1) * amp[0]; 
  jamp[4] += -thrust::complex<double> (0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 12
  VVV1_0(w[5], w[1], w[6], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[0]);
  jamp[1] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[3] += +amp[0]; 
  jamp[4] += -amp[0]; 
  // Amplitude(s) for diagram number 13
  FFV1_0(w[8], w[11], w[0], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[2] += -amp[0]; 
  // Amplitude(s) for diagram number 14
  FFV1_0(w[9], w[7], w[0], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[4] += -amp[0]; 
  // Amplitude(s) for diagram number 15
  VVV1_0(w[0], w[10], w[6], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[0]);
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[5] += +amp[0]; 
  VVVV1P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[10]);
  VVVV3P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[6]);
  VVVV4P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[9]);
  // Amplitude(s) for diagram number 16
  FFV1_0(w[3], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[5] += +amp[0]; 
  FFV1_0(w[3], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[1] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[4] += +amp[0]; 
  FFV1_0(w[3], w[2], w[9], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  jamp[0] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[4] += +amp[0]; 
  jamp[5] += -amp[0]; 
  // double CPPProcess::matrix_1_gg_ttxg() {


  int i, j; 
  // Local variables

  // const int ngraphs = 2;

  thrust::complex<double> ztemp; 
  // The color matrix;
  static const double denom[ncolor] = {9, 9, 9, 9, 9, 9}; 
  static const double cf[ncolor][ncolor] = {{64, -8, -8, 1, 1, 10}, {-8, 64, 1,
      10, -8, 1}, {-8, 1, 64, -8, 10, 1}, {1, 10, -8, 64, 1, -8}, {1, -8, 10,
      1, 64, -8}, {10, 1, 1, -8, -8, 64}};


  // Sum and square the color flows to get the matrix element
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + (ztemp * conj(jamp[i])).real()/denom[i]; 
  }

  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

}



CPPProcess::CPPProcess(int numiterations, int gpublocks, int gputhreads, 
bool verbose, bool debug)
: m_numiterations(numiterations), gpu_nblocks(gpublocks), 
gpu_nthreads(gputhreads), m_verbose(verbose), m_debug(debug), 
dim(gpu_nblocks * gpu_nthreads) 
{


  // Helicities for the process - nodim
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1, -1, -1,
      -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1, -1},
      {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1, -1,
      -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1}, {-1,
      1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1, -1}, {1, -1, -1,
      1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1, 1, -1}, {1, -1,
      1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1, -1, 1, -1}, {1,
      1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1, 1, 1, -1}, {1,
      1, 1, 1, 1}};
  gpuErrchk3(cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int)); )
  // perm - nodim
  // static int perm[nexternal] = {0, 1, 2, 3};
}

CPPProcess::~CPPProcess() {}

const std::vector<double> &CPPProcess::getMasses() const {return mME;}

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance(); 
  SLHAReader slha(param_card_name, m_verbose); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  if (m_verbose)
  {
    pars->printIndependentParameters(); 
    pars->printIndependentCouplings(); 
  }
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->ZERO); 
  static thrust::complex<double> tIPC[3] = {pars->GC_10, pars->GC_11,
      pars->GC_12};
  static double tIPD[2] = {pars->mdl_MT, pars->mdl_WT}; 

  gpuErrchk3(cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(thrust::complex<double>
      )));
  gpuErrchk3(cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(double))); 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

__global__ void sigmaKin(double * allmomenta, double * output) 
{
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();

  // Reset color flows

  // for (int xx = 0; xx < 384; ++xx) {
  const int nprocesses = 1; 
  int tid = blockIdx.x * blockDim.x + threadIdx.x; 

  // char *devPtr = (char *)tp.ptr;
  // size_t dpt = tp.pitch;
  // size_t slicePitch = dpt * 5;

  // char *dps = devPtr + dim * slicePitch;
  double matrix_element[nprocesses]; 

  thrust::complex<double> amp[1];  // was 18];


  double local_m[5][3]; 
  int DIM = blockDim.x * gridDim.x; 
  // for (int i=0; i<20;i++){
  // printf(" %f ", allmomenta[i]);
  // }
  // printf("\n");
  // printf("DIM is %i/%i\n", tid, DIM);
  for (int i = 0; i < 5; i++ )
  {
    for (int j = 0; j < 3; j++ )
    {
      local_m[i][j] = allmomenta[i * 3 * DIM + j * DIM + tid]; 
      // printf(" %f ", local_m[i][j]);
    }
    // printf("\n");
  }


  // Local variables and constants
  const int ncomb = 32; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  // double t[1];
  // Helicities for the process
  // static const int helicities[ncomb][nexternal] =
  // {{-1,-1,-1,-1,-1},{-1,-1,-1,-1,1},{-1,-1,-1,1,-1},{-1,-1,-1,1,1},{-1,-1,1,-
  // 1,-1},{-1,-1,1,-1,1},{-1,-1,1,1,-1},{-1,-1,1,1,1},{-1,1,-1,-1,-1},{-1,1,-1,
  // -1,1},{-1,1,-1,1,-1},{-1,1,-1,1,1},{-1,1,1,-1,-1},{-1,1,1,-1,1},{-1,1,1,1,-
  // 1},{-1,1,1,1,1},{1,-1,-1,-1,-1},{1,-1,-1,-1,1},{1,-1,-1,1,-1},{1,-1,-1,1,1}
  // ,{1,-1,1,-1,-1},{1,-1,1,-1,1},{1,-1,1,1,-1},{1,-1,1,1,1},{1,1,-1,-1,-1},{1,
  // 1,-1,-1,1},{1,1,-1,1,-1},{1,1,-1,1,1},{1,1,1,-1,-1},{1,1,1,-1,1},{1,1,1,1,-
  // 1},{1,1,1,1,1}};
  // Denominators: spins, colors and identical particles
  const int denominators[1] = {256}; 


  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  // int perm[nexternal];
  // for(int i = 0; i < nexternal; i++){
  // perm[i]=i;
  // }


  for (int ihel = 0; ihel < ncomb; ihel++ )
  {
    calculate_wavefunctions(ihel, local_m, matrix_element[0]); 
  }


  for (int i = 0; i < nprocesses; ++ i)
  {
    matrix_element[i] /= denominators[i]; 
  }
  for (int i = 0; i < nprocesses; ++ i)
  {
    output[i * nprocesses + tid] = matrix_element[i]; 
    // printf("output %i %i %i %f", tid, i, i*nprocesses+tid,
    // output[i*nprocesses+tid]);

  }


}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


