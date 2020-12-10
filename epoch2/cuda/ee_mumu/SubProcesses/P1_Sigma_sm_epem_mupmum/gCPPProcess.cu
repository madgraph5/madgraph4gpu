//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.8.2, 2020-10-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================


#include "../../src/HelAmps_sm.cu"


#include <algorithm>
#include <iostream>
#include "mgOnGpuTypes.h"
#include "mgOnGpuConfig.h"

#include "gCPPProcess.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > mu+ mu- WEIGHTED<=4 @1

#ifdef __CUDACC__
namespace gProc
#else
namespace Proc
#endif
{

using mgOnGpu::np4;  // 4: the dimension of 4-momenta (E,px,py,pz)
using mgOnGpu::npar;  // number of particles in total (initial + final)
using mgOnGpu::ncomb;  // number of helicity combinations

#ifdef __CUDACC__
__device__ __constant__ int cHel[ncomb][npar]; 
__device__ __constant__ fptype cIPC[6]; 
__device__ __constant__ fptype cIPD[2]; 
__device__ __constant__ int cNGoodHel[1]; 
__device__ __constant__ int cGoodHel[ncomb]; 
#else
static int cHel[ncomb][npar]; 
static fptype cIPC[6]; 
static fptype cIPD[2]; 
#endif

//--------------------------------------------------------------------------
using mgOnGpu::nwf; 
using mgOnGpu::nw6; 
//--------------------------------------------------------------------------

// Evaluate |M|^2 for each subprocess
// NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running sum
// of |M|^2 over helicities for the given event

__device__ void calculate_wavefunctions(int ihel, const fptype * allmomenta,
    fptype &meHelSum
#ifndef __CUDACC__
, const int ievt
#endif
)
{
  using namespace MG5_sm; 
  mgDebug(0, __FUNCTION__); 
  cxtype amp[1];  // was 2
  const int ncolor = 1; 
  cxtype jamp[ncolor]; 
  // Calculate wavefunctions for all processes
  using namespace MG5_sm; 
  cxtype w[nwf][nw6]; 
  for(int i = 0; i < 1; i++ )
  {
    jamp[i] = cxtype(0., 0.); 
  }

#ifdef __CUDACC__
  opzxxx(allmomenta, cHel[ihel][0], -1, w[0], 0); 
#else
  opzxxx(allmomenta, cHel[ihel][0], -1, w[0], ievt, 0); 
#endif 


#ifdef __CUDACC__
  imzxxx(allmomenta, cHel[ihel][1], +1, w[1], 1); 
#else
  imzxxx(allmomenta, cHel[ihel][1], +1, w[1], ievt, 1); 
#endif 


#ifdef __CUDACC__
  ixzxxx(allmomenta, cHel[ihel][2], -1, w[2], 2); 
#else
  ixzxxx(allmomenta, cHel[ihel][2], -1, w[2], ievt, 2); 
#endif 


#ifdef __CUDACC__
  oxzxxx(allmomenta, cHel[ihel][3], +1, w[3], 3); 
#else
  oxzxxx(allmomenta, cHel[ihel][3], +1, w[3], ievt, 3); 
#endif 

  FFV1P0_3(w[1], w[0], cxtype(cIPC[0], cIPC[1]), 0., 0., w[4]); 
  // Amplitude(s) for diagram number 1
  FFV1_0(w[2], w[3], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += -amp[0]; 
  FFV2_4_3(w[1], w[0], cxtype(cIPC[2], cIPC[3]), cxtype(cIPC[4], cIPC[5]),
      cIPD[0], cIPD[1], w[4]);
  // Amplitude(s) for diagram number 2
  FFV2_4_0(w[2], w[3], w[4], cxtype(cIPC[2], cIPC[3]), cxtype(cIPC[4],
      cIPC[5]), &amp[0]);
  jamp[0] += -amp[0]; 
  // double CPPProcess::matrix_1_epem_mupmum() {



  // Local variables

  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{1}}; 


  // Sum and square the color flows to get the matrix element
  for(int icol = 0; icol < ncolor; icol++ )
  {
    cxtype ztemp = cxmake(0, 0); 
    for(int jcol = 0; jcol < ncolor; jcol++ )
      ztemp = ztemp + cf[icol][jcol] * jamp[jcol]; 
    meHelSum = meHelSum + cxreal(ztemp * conj(jamp[icol]))/denom[icol]; 
  }

  // Store the leading color flows for choice of color
  // for(i=0;i < ncolor; i++)
  // jamp2[0][i] += real(jamp[i]*conj(jamp[i]));

  mgDebug(1, __FUNCTION__); 
  return; 
}



CPPProcess::CPPProcess(int numiterations, int gpublocks, int gputhreads, 
bool verbose, bool debug)
: m_numiterations(numiterations), gpu_nblocks(gpublocks), 
gpu_nthreads(gputhreads), m_verbose(verbose), 
dim(gpu_nblocks * gpu_nthreads) 
{


  // Helicities for the process - nodim
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1, -1, 1},
      {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1, -1, 1}, {-1, 1,
      1, -1}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, -1, 1, -1},
      {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 1}, {1, 1, 1, -1}, {1, 1, 1,
      1}};
#ifdef __CUDACC__
  checkCuda(cudaMemcpyToSymbol(cHel, tHel, ncomb * nexternal * sizeof(int))); 
#else
  memcpy(cHel, tHel, ncomb * nexternal * sizeof(int)); 
#endif

  // SANITY CHECK: GPU memory usage may be based on casts of fptype[2] to cxtype
  assert(sizeof(cxtype) == 2 * sizeof(fptype)); 
}

CPPProcess::~CPPProcess() {}

const std::vector<fptype> &CPPProcess::getMasses() const {return mME;}

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
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 

  static cxtype tIPC[3] = {pars->GC_3, pars->GC_50, pars->GC_59}; 
  static double tIPD[2] = {pars->mdl_MZ, pars->mdl_WZ}; 

#ifdef __CUDACC__
  checkCuda(cudaMemcpyToSymbol(cIPC, tIPC, 3 * sizeof(cxtype))); 
  checkCuda(cudaMemcpyToSymbol(cIPD, tIPD, 2 * sizeof(fptype))); 
#else
  memcpy(cIPC, tIPC, 3 * sizeof(cxtype)); 
  memcpy(cIPD, tIPD, 2 * sizeof(fptype)); 
#endif

}

//--------------------------------------------------------------------------

#ifdef __CUDACC__
__global__
void sigmaKin_getGoodHel(const fptype * allmomenta,  // input: momenta as AOSOA[npagM][npar][4][neppM] with nevt=npagM*neppM
bool * isGoodHel)  // output: isGoodHel[ncomb] - device array
{
  const int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
  fptype meHelSum[nprocesses] = {0};  // all zeros
  fptype meHelSumLast = 0; 
  for (int ihel = 0; ihel < ncomb; ihel++ )
  {
    // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running
    // sum of |M|^2 over helicities for the given event
    calculate_wavefunctions(ihel, allmomenta, meHelSum[0]); 
    if (meHelSum[0] != meHelSumLast)
    {
      isGoodHel[ihel] = true; 
      meHelSumLast = meHelSum[0]; 
    }
  }
}
#endif

//--------------------------------------------------------------------------

#ifdef __CUDACC__
void sigmaKin_setGoodHel(const bool * isGoodHel)  // input: isGoodHel[ncomb] - host array
{
  int nGoodHel[1] = {0}; 
  int goodHel[ncomb] = {0}; 
  for (int ihel = 0; ihel < ncomb; ihel++ )
  {
    // std::cout << "sigmaKin_setGoodHel ihel=" << ihel << ( isGoodHel[ihel] ?
    // " true" : " false" ) << std::endl;
    if (isGoodHel[ihel])
    {
      goodHel[nGoodHel[0]] = ihel; 
      nGoodHel[0]++; 
    }
  }
  checkCuda(cudaMemcpyToSymbol(cNGoodHel, nGoodHel, sizeof(int))); 
  checkCuda(cudaMemcpyToSymbol(cGoodHel, goodHel, ncomb * sizeof(int))); 
}
#endif



//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

__global__ void sigmaKin(const fptype * allmomenta, fptype * allMEs
#ifndef __CUDACC__
, const int nevt  // input: #events (for cuda: nevt == ndim == gpublocks*gputhreads)
#endif
) 
{
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();

#ifndef __CUDACC__
  const int maxtry = 10; 
  static unsigned long long sigmakin_itry = 0;  // first iteration over nevt events
  static bool sigmakin_goodhel[ncomb] = {false}; 
#endif

  // Reset color flows


  // start sigmakin_lines


  mgDebugInitialise(); 
  // Set the parameters which change event by event
  // Need to discuss this with Stefan
  // pars->setDependentParameters();
  // pars->setDependentCouplings();
  // Reset color flows

#ifndef __CUDACC__
  //** START LOOP ON IEVT **
  for (int ievt = 0; ievt < nevt; ++ ievt)
#endif
  {
#ifdef __CUDACC__
    const int idim = blockDim.x * blockIdx.x + threadIdx.x;  // event# == threadid (previously was: tid)
    const int ievt = idim; 
    // printf( "sigmakin: ievt %d\n", ievt );
#endif

    // Denominators: spins, colors and identical particles
    const int nprocesses = 1;  // FIXME: assume process.nprocesses == 1
    const int denominators[1] = {4}; 

    // Reset the "matrix elements" - running sums of |M|^2 over helicities for
    // the given event
    fptype meHelSum[nprocesses] = {0};  // all zeros

#ifdef __CUDACC__
    // CUDA - using precomputed good helicities
    for (int ighel = 0; ighel < cNGoodHel[0]; ighel++ )
    {
      const int ihel = cGoodHel[ighel]; 
      calculate_wavefunctions(ihel, allmomenta, meHelSum[0]); 
    }
#else
    // C++ - compute good helicities within this loop
    fptype meHelSumLast = 0;  // check for good helicities
    for (int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (sigmakin_itry > maxtry && !sigmakin_goodhel[ihel])
        continue; 
      // NB: calculate_wavefunctions ADDS |M|^2 for a given ihel to the running
      // sum of |M|^2 over helicities for the given event
      calculate_wavefunctions(ihel, allmomenta, meHelSum[0], ievt); 
      if (sigmakin_itry <= maxtry)
      {
        if ( !sigmakin_goodhel[ihel] && meHelSum[0] > meHelSumLast)
          sigmakin_goodhel[ihel] = true; 
        meHelSumLast = meHelSum[0]; 
      }
    }
#endif

    // Get the final |M|^2 as an average over helicities/colors of the running
    // sum of |M|^2 over helicities for the given event
    // [NB 'sum over final spins, average over initial spins', eg see
    // https://www.uzh.ch/cmsssl/physik/dam/jcr:2e24b7b1-f4d7-4160-817e-47b13dbf
    // 1d7c/Handout_4_2016-UZH.pdf]
    for (int iproc = 0; iproc < nprocesses; ++ iproc)
    {
      meHelSum[iproc] /= denominators[iproc]; 
    }

    // Set the final average |M|^2 for this event in the output array for all
    // events
    for (int iproc = 0; iproc < nprocesses; ++ iproc)
    {
      allMEs[iproc * nprocesses + ievt] = meHelSum[iproc]; 
    }

#ifndef __CUDACC__
    if (sigmakin_itry <= maxtry)
      sigmakin_itry++; 
    // if ( sigmakin_itry == maxtry )
    // for (int ihel = 0; ihel < ncomb; ihel++ )
    // printf( "sigmakin: ihelgood %2d %d\n", ihel, sigmakin_goodhel[ihel] );
#endif
    //** END LOOP ON IEVT **
    mgDebugFinalise(); 

  }

  //--------------------------------------------------------------------------

}



}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------


