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
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

__constant__ int cHel[16][4]; 
// __constant__ double cmME[4]; value hardcoded now
// extern __constant__ int cPerm[4];
// 
__constant__ double cIPC[6];  // coupling ?
__constant__ double cIPD[2]; 


__device__ void printComplex6( const thrust::complex<double> c6[6], const char* txt )
{
  for (int i=0; i<6; i++)
    printf("%s[%i]=(%f,%f)\n", txt, i, c6[i].real(), c6[i].imag());
}

__device__ void printDouble4( const double p[4], const char* txt )
{
  //for (int i=0; i<4; i++) printf("%s[%i]=%f\n", txt, i, p[i]);
  printf("%s[0..3] = ", txt);
  for (int i=0; i<4; i++) printf("%f ", p[i]);
  printf("\n");
}

// Evaluate |M|^2 for each subprocess

__device__ void calculate_wavefunctions(int ihel, double local_mom[4][3],
    double &matrix)
{
  thrust::complex<double> amp[2]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[5][6]; 
  oxxxxx(local_mom[0], 0., cHel[ihel][0], -1, w[0]); 
  ixxxxx(local_mom[1], 0., cHel[ihel][1], +1, w[1]); 
  ixxxxx(local_mom[2], 0., cHel[ihel][2], -1, w[2]); 
  oxxxxx(local_mom[3], 0., cHel[ihel][3], +1, w[3]); 
  FFV1P0_3(w[1], w[0], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[4]);
  // Amplitude(s) for diagram number 1
  FFV1_0(w[2], w[3], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[0]);
  FFV2_4_3(w[1], w[0], thrust::complex<double> (cIPC[2], cIPC[3]),
      thrust::complex<double> (cIPC[4], cIPC[5]), cIPD[0], cIPD[1], w[4]);
  // Amplitude(s) for diagram number 2
  FFV2_4_0(w[2], w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      thrust::complex<double> (cIPC[4], cIPC[5]), &amp[1]);

  //printf("DEBUG0: amp0=(%f,%f) amp1=(%f,%f)\n", amp[0].real(), amp[0].imag(), amp[1].real(), amp[1].imag());

  // double CPPProcess::matrix_1_epem_mupmum() {
  int i, j; 
  // Local variables

  // const int ngraphs = 2;
  const int ncolor = 1; 
  thrust::complex<double> ztemp; 
  thrust::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{1}}; 

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1]; 

  // Sum and square the color flows to get the matrix element
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + (ztemp * conj(jamp[i])).real()/denom[i]; 
  }

  //printf("DEBUG0: ME=%f\n", matrix);

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
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  static thrust::complex<double> tIPC[3] = {pars->GC_3, pars->GC_50,
      pars->GC_59};
  static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ}; 

  cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(thrust::complex<double> )); 
  cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(double)); 
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
  // size_t slicePitch = dpt * 4;

  // char *dps = devPtr + dim * slicePitch;
  double matrix_element[nprocesses]; 

  thrust::complex<double> amp[2]; 

  double local_m[4][3]; 
  int DIM = blockDim.x * gridDim.x; 
  // for (int i=0; i<20;i++){
  // printf(" %f ", allmomenta[i]);
  // }
  // printf("\n");
  // printf("DIM is %i/%i\n", tid, DIM);
  for (int i = 0; i < 4; i++ )
  {
    for (int j = 0; j < 3; j++ )
    {
      local_m[i][j] = allmomenta[i * 3 * DIM + j * DIM + tid]; 
      // printf(" %f ", local_m[i][j]);
    }
    // printf("\n");
  }


  // Local variables and constants
  const int ncomb = 16; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  // double t[1];
  // Helicities for the process
  // static const int helicities[ncomb][nexternal] =
  // {{-1,-1,-1,-1},{-1,-1,-1,1},{-1,-1,1,-1},{-1,-1,1,1},{-1,1,-1,-1},{-1,1,-1,
  // 1},{-1,1,1,-1},{-1,1,1,1},{1,-1,-1,-1},{1,-1,-1,1},{1,-1,1,-1},{1,-1,1,1},{
  // 1,1,-1,-1},{1,1,-1,1},{1,1,1,-1},{1,1,1,1}};
  // Denominators: spins, colors and identical particles
  const int denominators[1] = {4}; 


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


  //printf("\n\nblock (%i / %i), thread (%i)\n\n", blockIdx.x, blockDim.x, threadIdx.x);
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
  //printf("DEBUG: ME=%f\n", matrix_element[0]);

}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


