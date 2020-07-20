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
// Process: g g > t t~ g WEIGHTED<=3 @1

__constant__ int cHel[32][5]; 
// __constant__ double cmME[5]; value hardcoded now
// extern __constant__ int cPerm[4];
// 
__constant__ double cIPC[6];  // coupling ?
__constant__ double cIPD[2]; 


// Evaluate |M|^2 for each subprocess

__device__ void calculate_wavefunctions(int ihel, char * dps, size_t dpt,
    double &matrix)
{
  thrust::complex<double> amp[18]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[12][6]; 
  vxxxxx((double * )(dps + 0 * dpt), 0., cHel[ihel][0], -1, w[0]); 
  vxxxxx((double * )(dps + 1 * dpt), 0., cHel[ihel][1], -1, w[1]); 
  oxxxxx((double * )(dps + 2 * dpt), cIPD[0], cHel[ihel][2], +1, w[2]); 
  ixxxxx((double * )(dps + 3 * dpt), cIPD[0], cHel[ihel][3], -1, w[3]); 
  vxxxxx((double * )(dps + 4 * dpt), 0., cHel[ihel][4], +1, w[4]); 
  VVV1P0_1(w[0], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[5]);
  FFV1P0_3(w[3], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[6]);
  // Amplitude(s) for diagram number 1
  VVV1_0(w[5], w[6], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[0]); printf(" &amp[0] %f %f", &amp[0].real(), &amp[0].imag(); )
  FFV1_1(w[2], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[7]);
  // Amplitude(s) for diagram number 2
  FFV1_0(w[3], w[7], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[1]); printf(" &amp[1] %f %f", &amp[1].real(), &amp[1].imag(); )
  FFV1_2(w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[8]);
  // Amplitude(s) for diagram number 3
  FFV1_0(w[8], w[2], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[2]); printf(" &amp[2] %f %f", &amp[2].real(), &amp[2].imag(); )
  FFV1_1(w[2], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[5]);
  FFV1_2(w[3], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[9]);
  // Amplitude(s) for diagram number 4
  FFV1_0(w[9], w[5], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[3]); printf(" &amp[3] %f %f", &amp[3].real(), &amp[3].imag(); )
  VVV1P0_1(w[1], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[10]);
  // Amplitude(s) for diagram number 5
  FFV1_0(w[3], w[5], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[4]); printf(" &amp[4] %f %f", &amp[4].real(), &amp[4].imag(); )
  // Amplitude(s) for diagram number 6
  FFV1_0(w[8], w[5], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[5]); printf(" &amp[5] %f %f", &amp[5].real(), &amp[5].imag(); )
  FFV1_2(w[3], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[5]);
  FFV1_1(w[2], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[11]);
  // Amplitude(s) for diagram number 7
  FFV1_0(w[5], w[11], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[6]); printf(" &amp[6] %f %f", &amp[6].real(), &amp[6].imag(); )
  // Amplitude(s) for diagram number 8
  FFV1_0(w[5], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[7]); printf(" &amp[7] %f %f", &amp[7].real(), &amp[7].imag(); )
  // Amplitude(s) for diagram number 9
  FFV1_0(w[5], w[7], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[8]); printf(" &amp[8] %f %f", &amp[8].real(), &amp[8].imag(); )
  VVV1P0_1(w[0], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[5]);
  // Amplitude(s) for diagram number 10
  FFV1_0(w[3], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[9]); printf(" &amp[9] %f %f", &amp[9].real(), &amp[9].imag(); )
  // Amplitude(s) for diagram number 11
  FFV1_0(w[9], w[2], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[10]); printf(" &amp[10] %f %f", &amp[10].real(), &amp[10].imag(); )
  // Amplitude(s) for diagram number 12
  VVV1_0(w[5], w[1], w[6], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[11]); printf(" &amp[11] %f %f", &amp[11].real(), &amp[11].imag(); )
  // Amplitude(s) for diagram number 13
  FFV1_0(w[8], w[11], w[0], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[12]); printf(" &amp[12] %f %f", &amp[12].real(), &amp[12].imag(); )
  // Amplitude(s) for diagram number 14
  FFV1_0(w[9], w[7], w[0], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[13]); printf(" &amp[13] %f %f", &amp[13].real(), &amp[13].imag(); )
  // Amplitude(s) for diagram number 15
  VVV1_0(w[0], w[10], w[6], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[14]); printf(" &amp[14] %f %f", &amp[14].real(), &amp[14].imag(); )
  VVVV1P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[10]);
  VVVV3P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[6]);
  VVVV4P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[9]);
  // Amplitude(s) for diagram number 16
  FFV1_0(w[3], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[15]); printf(" &amp[15] %f %f", &amp[15].real(), &amp[15].imag(); )
  FFV1_0(w[3], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[16]); printf(" &amp[16] %f %f", &amp[16].real(), &amp[16].imag(); )
  FFV1_0(w[3], w[2], w[9], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[17]); printf(" &amp[17] %f %f", &amp[17].real(), &amp[17].imag(); )
  // double CPPProcess::matrix_1_gg_ttxg() {
  int i, j; 
  // Local variables

  // const int ngraphs = 2;
  const int ncolor = 6; 
  thrust::complex<double> ztemp; 
  thrust::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {9, 9, 9, 9, 9, 9}; 
  static const double cf[ncolor][ncolor] = {{64, -8, -8, 1, 1, 10}, {-8, 64, 1,
      10, -8, 1}, {-8, 1, 64, -8, 10, 1}, {1, 10, -8, 64, 1, -8}, {1, -8, 10,
      1, 64, -8}, {10, 1, 1, -8, -8, 64}};

  // Calculate color flows
  jamp[0] = -amp[0] + thrust::complex<double> (0, 1) * amp[2] +
      thrust::complex<double> (0, 1) * amp[4] - amp[5] + amp[14] + amp[15] -
      amp[17];
  jamp[1] = -amp[3] - thrust::complex<double> (0, 1) * amp[4] +
      thrust::complex<double> (0, 1) * amp[10] + amp[11] - amp[14] - amp[15] -
      amp[16];
  jamp[2] = +amp[0] - thrust::complex<double> (0, 1) * amp[2] +
      thrust::complex<double> (0, 1) * amp[9] - amp[11] - amp[12] + amp[16] +
      amp[17];
  jamp[3] = -amp[6] + thrust::complex<double> (0, 1) * amp[7] -
      thrust::complex<double> (0, 1) * amp[9] + amp[11] - amp[14] - amp[15] -
      amp[16];
  jamp[4] = +amp[0] + thrust::complex<double> (0, 1) * amp[1] -
      thrust::complex<double> (0, 1) * amp[10] - amp[11] - amp[13] + amp[16] +
      amp[17];
  jamp[5] = -amp[0] - thrust::complex<double> (0, 1) * amp[1] -
      thrust::complex<double> (0, 1) * amp[7] - amp[8] + amp[14] + amp[15] -
      amp[17];

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
  mME.push_back(pars->ZERO); 
  static thrust::complex<double> tIPC[3] = {pars->GC_10, pars->GC_11,
      pars->GC_12};
  static double tIPD[2] = {pars->mdl_MT, pars->mdl_WT}; 

  cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(thrust::complex<double> )); 
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

  thrust::complex<double> amp[18]; 


  // Local variables and constants
  const int ncomb = 32; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  double t[1]; 
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


