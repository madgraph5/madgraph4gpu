//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3.py3, 2020-06-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "HelAmps_sm.h"

#include <algorithm> 
#include <iostream> 
#include <thrust/complex.h> 

using namespace MG5_sm; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > t t~ WEIGHTED<=2 @1

__constant__ int cHel[16][4]; 
// __constant__ double cmME[4]; value hardcoded now
// extern __constant__ int cPerm[4];
// 
__constant__ double cIPC[4];  // coupling ?
__constant__ double cIPD[2]; 


// Evaluate |M|^2 for each subprocess

__device__ void calculate_wavefunctions(int ihel, char * dps, size_t dpt,
    double &matrix)
{
  thrust::complex<double> amp[3]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[5][6]; 
  vxxxxx((double * )(dps + 0 * dpt), 0., cHel[ihel][0], -1, w[0]); 
  vxxxxx((double * )(dps + 1 * dpt), 0., cHel[ihel][1], -1, w[1]); 
  oxxxxx((double * )(dps + 2 * dpt), cIPD[0], cHel[ihel][2], +1, w[2]); 
  ixxxxx((double * )(dps + 3 * dpt), cIPD[0], cHel[ihel][3], -1, w[3]); 
  VVV1P0_1(w[0], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[4]);
  // Amplitude(s) for diagram number 1
  FFV1_0(w[3], w[2], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[0]);
  FFV1_1(w[2], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[4]);
  // Amplitude(s) for diagram number 2
  FFV1_0(w[3], w[4], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[1]);
  FFV1_2(w[3], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[4]);
  // Amplitude(s) for diagram number 3
  FFV1_0(w[4], w[2], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[2]);
  // double CPPProcess::matrix_1_gg_ttx() {
  int i, j; 
  // Local variables

  // const int ngraphs = 2;
  const int ncolor = 2; 
  thrust::complex<double> ztemp; 
  thrust::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {3, 3}; 
  static const double cf[ncolor][ncolor] = {{16, -2}, {-2, 16}}; 

  // Calculate color flows
  jamp[0] = +thrust::complex<double> (0, 1) * amp[0] - amp[1]; 
  jamp[1] = -thrust::complex<double> (0, 1) * amp[0] - amp[2]; 

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
gpu_nthreads(gputhreads), dim(gpu_nblocks * gpu_nthreads) 
{


  // Helicities for the process - nodim
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1, -1, 1},
      {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1}, {-1, 1,
      1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1, 1, -1},
      {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1, 1, 1,
      1}};
  cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int)); 
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
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->mdl_MT); 
  static thrust::complex<double> tIPC[2] = {pars->GC_10, pars->GC_11}; 
  static double tIPD[2] = {pars->mdl_MT, pars->mdl_WT}; 

  cudaMemcpyToSymbol(cIPC, tIPC, 2 * sizeof(thrust::complex<double> )); 
  cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(double)); 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

__global__ void sigmaKin(cudaPitchedPtr tp, double * meDevPtr, size_t mePitch) 
{
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();

  // Reset color flows

  // for (int xx = 0; xx < 384; ++xx) {
  int nprocesses = 1; 
  int dim = blockIdx.x * blockDim.x + threadIdx.x; 

  char * devPtr = (char * )tp.ptr; 
  size_t dpt = tp.pitch; 
  size_t slicePitch = dpt * 4; 

  char * dps = devPtr + dim * slicePitch; 
  double * matrix_element = (double * )((char * )meDevPtr + dim * mePitch); 

  thrust::complex<double> amp[3]; 


  // Local variables and constants
  const int ncomb = 16; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  double t[1]; 
  // Helicities for the process
  // static const int helicities[ncomb][nexternal] =
  // {{-1,-1,-1,-1},{-1,-1,-1,1},{-1,-1,1,-1},{-1,-1,1,1},{-1,1,-1,-1},{-1,1,-1,
  // 1},{-1,1,1,-1},{-1,1,1,1},{1,-1,-1,-1},{1,-1,-1,1},{1,-1,1,-1},{1,-1,1,1},{
  // 1,1,-1,-1},{1,1,-1,1},{1,1,1,-1},{1,1,1,1}};
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
    calculate_wavefunctions(ihel, dps, dpt, matrix_element[0]); 
  }


  for (int i = 0; i < nprocesses; ++ i)
  {
    matrix_element[i] /= denominators[i]; 
  }


}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


