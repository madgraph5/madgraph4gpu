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
// Process: g g > t t~ g g WEIGHTED<=4 @1

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
  cxtype amp[1];  // was 159
  const int ncolor = 24; 
  cxtype jamp[ncolor]; 
  // Calculate wavefunctions for all processes
  using namespace MG5_sm; 
  cxtype w[nwf][nw6]; 
  for(int i = 0; i < 24; i++ )
  {
    jamp[i] = cxtype(0., 0.); 
  }

#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][0], -1, w[0], 0); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][0], -1, w[0], ievt, 0); 
#endif 


#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][1], -1, w[1], 1); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][1], -1, w[1], ievt, 1); 
#endif 


#ifdef __CUDACC__
  oxxxxx(allmomenta, cIPD[0], cHel[ihel][2], +1, w[2], 2); 
#else
  oxxxxx(allmomenta, cIPD[0], cHel[ihel][2], +1, w[2], ievt, 2); 
#endif 


#ifdef __CUDACC__
  ixxxxx(allmomenta, cIPD[0], cHel[ihel][3], -1, w[3], 3); 
#else
  ixxxxx(allmomenta, cIPD[0], cHel[ihel][3], -1, w[3], ievt, 3); 
#endif 


#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][4], +1, w[4], 4); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][4], +1, w[4], ievt, 4); 
#endif 


#ifdef __CUDACC__
  vxxxxx(allmomenta, 0., cHel[ihel][5], +1, w[5], 5); 
#else
  vxxxxx(allmomenta, 0., cHel[ihel][5], +1, w[5], ievt, 5); 
#endif 

  VVV1P0_1(w[0], w[1], cxtype(cIPC[0], cIPC[1]), 0., 0., w[6]); 
  FFV1P0_3(w[3], w[2], cxtype(cIPC[2], cIPC[3]), 0., 0., w[7]); 
  // Amplitude(s) for diagram number 1
  VVVV1_0(w[6], w[7], w[4], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[6], w[7], w[4], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[6], w[7], w[4], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[6], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[8]); 
  // Amplitude(s) for diagram number 2
  VVV1_0(w[7], w[5], w[8], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[6], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[9]); 
  // Amplitude(s) for diagram number 3
  VVV1_0(w[7], w[4], w[9], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[4], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[10]); 
  // Amplitude(s) for diagram number 4
  VVV1_0(w[6], w[7], w[10], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[11]); 
  FFV1_2(w[3], w[6], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  // Amplitude(s) for diagram number 5
  FFV1_0(w[12], w[11], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 6
  FFV1_0(w[3], w[11], w[9], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +amp[0]; 
  jamp[14] += -amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[17] += +amp[0]; 
  FFV1_2(w[3], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[13]); 
  // Amplitude(s) for diagram number 7
  FFV1_0(w[13], w[11], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[14]); 
  // Amplitude(s) for diagram number 8
  FFV1_0(w[12], w[14], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 9
  FFV1_0(w[3], w[14], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +amp[0]; 
  jamp[20] += -amp[0]; 
  jamp[22] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_2(w[3], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[15]); 
  // Amplitude(s) for diagram number 10
  FFV1_0(w[15], w[14], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[6], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[16]); 
  // Amplitude(s) for diagram number 11
  FFV1_0(w[15], w[16], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 12
  FFV1_0(w[15], w[2], w[9], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[18] += -amp[0]; 
  jamp[20] += +amp[0]; 
  // Amplitude(s) for diagram number 13
  FFV1_0(w[13], w[16], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 14
  FFV1_0(w[13], w[2], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[6] += -amp[0]; 
  jamp[12] += -amp[0]; 
  jamp[14] += +amp[0]; 
  // Amplitude(s) for diagram number 15
  FFV1_0(w[3], w[16], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[6] += -amp[0]; 
  jamp[7] += +amp[0]; 
  // Amplitude(s) for diagram number 16
  FFV1_0(w[12], w[2], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += +amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[22] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_1(w[2], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  FFV1_2(w[3], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[16]); 
  FFV1_1(w[12], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[8]); 
  // Amplitude(s) for diagram number 17
  FFV1_0(w[16], w[8], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += -amp[0]; 
  FFV1_1(w[12], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[9]); 
  // Amplitude(s) for diagram number 18
  FFV1_0(w[16], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += -amp[0]; 
  // Amplitude(s) for diagram number 19
  FFV1_0(w[16], w[12], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[1], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[6]); 
  FFV1P0_3(w[3], w[12], cxtype(cIPC[2], cIPC[3]), 0., 0., w[17]); 
  // Amplitude(s) for diagram number 20
  VVV1_0(w[6], w[5], w[17], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[4] += -amp[0]; 
  jamp[5] += +amp[0]; 
  // Amplitude(s) for diagram number 21
  FFV1_0(w[3], w[9], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 22
  FFV1_0(w[13], w[12], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[1], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[18]); 
  // Amplitude(s) for diagram number 23
  VVV1_0(w[18], w[4], w[17], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[3] += +amp[0]; 
  jamp[4] += -amp[0]; 
  // Amplitude(s) for diagram number 24
  FFV1_0(w[3], w[8], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 25
  FFV1_0(w[15], w[12], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[12], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[19]); 
  // Amplitude(s) for diagram number 26
  FFV1_0(w[15], w[19], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += -amp[0]; 
  // Amplitude(s) for diagram number 27
  FFV1_0(w[15], w[9], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += -amp[0]; 
  // Amplitude(s) for diagram number 28
  FFV1_0(w[13], w[19], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += -amp[0]; 
  // Amplitude(s) for diagram number 29
  FFV1_0(w[13], w[8], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += -amp[0]; 
  // Amplitude(s) for diagram number 30
  FFV1_0(w[3], w[19], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 31
  VVV1_0(w[1], w[10], w[17], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[5] += +amp[0]; 
  VVVV1P0_1(w[1], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[17]); 
  VVVV3P0_1(w[1], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[19]); 
  VVVV4P0_1(w[1], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[8]); 
  // Amplitude(s) for diagram number 32
  FFV1_0(w[3], w[12], w[17], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[1] += -amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[5] += +amp[0]; 
  FFV1_0(w[3], w[12], w[19], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[4] += +amp[0]; 
  FFV1_0(w[3], w[12], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += -amp[0]; 
  jamp[2] += +amp[0]; 
  jamp[4] += +amp[0]; 
  jamp[5] += -amp[0]; 
  FFV1_2(w[3], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  FFV1_1(w[2], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[9]); 
  FFV1_2(w[12], w[4], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[20]); 
  // Amplitude(s) for diagram number 33
  FFV1_0(w[20], w[9], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[11] += -amp[0]; 
  FFV1_2(w[12], w[5], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[21]); 
  // Amplitude(s) for diagram number 34
  FFV1_0(w[21], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += -amp[0]; 
  // Amplitude(s) for diagram number 35
  FFV1_0(w[12], w[9], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  FFV1P0_3(w[12], w[2], cxtype(cIPC[2], cIPC[3]), 0., 0., w[22]); 
  // Amplitude(s) for diagram number 36
  VVV1_0(w[6], w[5], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[9] += +amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[23] += +amp[0]; 
  // Amplitude(s) for diagram number 37
  FFV1_0(w[21], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 38
  FFV1_0(w[12], w[14], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 39
  VVV1_0(w[18], w[4], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[11] += +amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[17] += +amp[0]; 
  jamp[21] += -amp[0]; 
  // Amplitude(s) for diagram number 40
  FFV1_0(w[20], w[2], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 41
  FFV1_0(w[12], w[11], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  FFV1_2(w[12], w[1], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[23]); 
  // Amplitude(s) for diagram number 42
  FFV1_0(w[23], w[11], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[17] += -amp[0]; 
  // Amplitude(s) for diagram number 43
  FFV1_0(w[21], w[11], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[15] += -amp[0]; 
  // Amplitude(s) for diagram number 44
  FFV1_0(w[23], w[14], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[23] += -amp[0]; 
  // Amplitude(s) for diagram number 45
  FFV1_0(w[20], w[14], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[21] += -amp[0]; 
  // Amplitude(s) for diagram number 46
  FFV1_0(w[23], w[2], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 47
  VVV1_0(w[1], w[10], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[9] += +amp[0]; 
  jamp[11] += -amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[23] += +amp[0]; 
  // Amplitude(s) for diagram number 48
  FFV1_0(w[12], w[2], w[17], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += +amp[0]; 
  jamp[11] += -amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_0(w[12], w[2], w[19], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[11] += -amp[0]; 
  jamp[15] += +amp[0]; 
  jamp[17] += -amp[0]; 
  jamp[21] += +amp[0]; 
  FFV1_0(w[12], w[2], w[8], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[9] += -amp[0]; 
  jamp[15] += +amp[0]; 
  jamp[21] += +amp[0]; 
  jamp[23] += -amp[0]; 
  VVV1P0_1(w[0], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[12]); 
  FFV1_2(w[3], w[12], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[22]); 
  // Amplitude(s) for diagram number 49
  FFV1_0(w[22], w[9], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[12], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[23]); 
  // Amplitude(s) for diagram number 50
  FFV1_0(w[3], w[9], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[11] += +amp[0]; 
  // Amplitude(s) for diagram number 51
  FFV1_0(w[13], w[9], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[12], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[20]); 
  // Amplitude(s) for diagram number 52
  FFV1_0(w[16], w[20], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 53
  FFV1_0(w[16], w[2], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[19] += -amp[0]; 
  jamp[22] += +amp[0]; 
  // Amplitude(s) for diagram number 54
  FFV1_0(w[16], w[14], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 55
  FFV1_0(w[3], w[20], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +amp[0]; 
  jamp[3] += -amp[0]; 
  jamp[12] += -amp[0]; 
  jamp[13] += +amp[0]; 
  // Amplitude(s) for diagram number 56
  FFV1_0(w[22], w[2], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += +amp[0]; 
  jamp[11] += -amp[0]; 
  jamp[20] += -amp[0]; 
  jamp[21] += +amp[0]; 
  // Amplitude(s) for diagram number 57
  VVV1_0(w[12], w[18], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 58
  VVVV1_0(w[12], w[1], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[12], w[1], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[12], w[1], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[12], w[1], cxtype(cIPC[0], cIPC[1]), 0., 0., w[21]); 
  // Amplitude(s) for diagram number 59
  VVV1_0(w[7], w[5], w[21], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 60
  VVV1_0(w[1], w[7], w[23], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 61
  FFV1_0(w[3], w[14], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += +amp[0]; 
  jamp[20] += -amp[0]; 
  jamp[21] += +amp[0]; 
  jamp[22] += -amp[0]; 
  // Amplitude(s) for diagram number 62
  FFV1_0(w[22], w[14], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 63
  FFV1_0(w[13], w[2], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +amp[0]; 
  jamp[6] += -amp[0]; 
  jamp[8] += +amp[0]; 
  jamp[12] += -amp[0]; 
  // Amplitude(s) for diagram number 64
  FFV1_0(w[13], w[20], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[5], cxtype(cIPC[0], cIPC[1]), 0., 0., w[20]); 
  FFV1_2(w[3], w[20], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[21]); 
  // Amplitude(s) for diagram number 65
  FFV1_0(w[21], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[20], w[4], cxtype(cIPC[0], cIPC[1]), 0., 0., w[22]); 
  // Amplitude(s) for diagram number 66
  FFV1_0(w[3], w[9], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += +amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[9] += +amp[0]; 
  jamp[10] += -amp[0]; 
  // Amplitude(s) for diagram number 67
  FFV1_0(w[15], w[9], w[20], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[2], w[20], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[23]); 
  // Amplitude(s) for diagram number 68
  FFV1_0(w[16], w[23], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 69
  FFV1_0(w[16], w[2], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[16] += +amp[0]; 
  jamp[19] += -amp[0]; 
  // Amplitude(s) for diagram number 70
  FFV1_0(w[16], w[11], w[20], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 71
  FFV1_0(w[3], w[23], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +amp[0]; 
  jamp[5] += -amp[0]; 
  jamp[18] += -amp[0]; 
  jamp[19] += +amp[0]; 
  // Amplitude(s) for diagram number 72
  FFV1_0(w[21], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += +amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[14] += -amp[0]; 
  jamp[15] += +amp[0]; 
  // Amplitude(s) for diagram number 73
  VVV1_0(w[20], w[6], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 74
  VVVV1_0(w[20], w[1], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[20], w[1], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[20], w[1], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[20], w[1], cxtype(cIPC[0], cIPC[1]), 0., 0., w[12]); 
  // Amplitude(s) for diagram number 75
  VVV1_0(w[7], w[4], w[12], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 76
  VVV1_0(w[1], w[7], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 77
  FFV1_0(w[3], w[11], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += +amp[0]; 
  jamp[14] += -amp[0]; 
  jamp[15] += +amp[0]; 
  jamp[16] += -amp[0]; 
  // Amplitude(s) for diagram number 78
  FFV1_0(w[21], w[11], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 79
  FFV1_0(w[15], w[2], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[10] += +amp[0]; 
  jamp[18] += -amp[0]; 
  // Amplitude(s) for diagram number 80
  FFV1_0(w[15], w[23], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  FFV1_1(w[9], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[23]); 
  // Amplitude(s) for diagram number 81
  FFV1_0(w[15], w[23], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += -amp[0]; 
  FFV1_2(w[15], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[12]); 
  // Amplitude(s) for diagram number 82
  FFV1_0(w[12], w[9], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += -amp[0]; 
  // Amplitude(s) for diagram number 83
  FFV1_0(w[13], w[23], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += -amp[0]; 
  FFV1_2(w[13], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[21]); 
  // Amplitude(s) for diagram number 84
  FFV1_0(w[21], w[9], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += -amp[0]; 
  // Amplitude(s) for diagram number 85
  FFV1_0(w[3], w[23], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[10], cxtype(cIPC[0], cIPC[1]), 0., 0., w[23]); 
  // Amplitude(s) for diagram number 86
  FFV1_0(w[3], w[9], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[11] += +amp[0]; 
  FFV1_2(w[16], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[22]); 
  // Amplitude(s) for diagram number 87
  FFV1_0(w[22], w[11], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += -amp[0]; 
  FFV1_1(w[11], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[20]); 
  // Amplitude(s) for diagram number 88
  FFV1_0(w[16], w[20], w[5], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += -amp[0]; 
  // Amplitude(s) for diagram number 89
  FFV1_0(w[22], w[14], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[22] += -amp[0]; 
  FFV1_1(w[14], w[0], cxtype(cIPC[2], cIPC[3]), cIPD[0], cIPD[1], w[24]); 
  // Amplitude(s) for diagram number 90
  FFV1_0(w[16], w[24], w[4], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += -amp[0]; 
  // Amplitude(s) for diagram number 91
  FFV1_0(w[22], w[2], w[10], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 92
  FFV1_0(w[16], w[2], w[23], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +amp[0]; 
  jamp[5] += -amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[22] += +amp[0]; 
  // Amplitude(s) for diagram number 93
  VVVV1_0(w[0], w[6], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[0], w[6], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[0], w[6], w[7], w[5], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[6], cxtype(cIPC[0], cIPC[1]), 0., 0., w[22]); 
  // Amplitude(s) for diagram number 94
  VVV1_0(w[7], w[5], w[22], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[7], cxtype(cIPC[0], cIPC[1]), 0., 0., w[25]); 
  // Amplitude(s) for diagram number 95
  VVV1_0(w[6], w[5], w[25], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 96
  FFV1_0(w[3], w[14], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +amp[0]; 
  jamp[19] += -amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[23] += +amp[0]; 
  // Amplitude(s) for diagram number 97
  FFV1_0(w[3], w[24], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 98
  FFV1_0(w[13], w[2], w[22], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[14] += +amp[0]; 
  // Amplitude(s) for diagram number 99
  FFV1_0(w[21], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 100
  VVVV1_0(w[0], w[18], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[0], w[18], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[0], w[18], w[7], w[4], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  VVV1P0_1(w[0], w[18], cxtype(cIPC[0], cIPC[1]), 0., 0., w[6]); 
  // Amplitude(s) for diagram number 101
  VVV1_0(w[7], w[4], w[6], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 102
  VVV1_0(w[18], w[4], w[25], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 103
  FFV1_0(w[3], w[11], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[17] += +amp[0]; 
  // Amplitude(s) for diagram number 104
  FFV1_0(w[3], w[20], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 105
  FFV1_0(w[15], w[2], w[6], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[4] += -amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[20] += +amp[0]; 
  // Amplitude(s) for diagram number 106
  FFV1_0(w[12], w[2], w[18], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 107
  VVVV1_0(w[0], w[1], w[7], w[10], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV3_0(w[0], w[1], w[7], w[10], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVVV4_0(w[0], w[1], w[7], w[10], cxtype(cIPC[4], cIPC[5]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 108
  VVV1_0(w[1], w[10], w[25], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 109
  VVV1_0(w[1], w[7], w[23], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 110
  FFV1_0(w[13], w[20], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += -amp[0]; 
  // Amplitude(s) for diagram number 111
  FFV1_0(w[21], w[11], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[14] += -amp[0]; 
  // Amplitude(s) for diagram number 112
  FFV1_0(w[15], w[24], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += -amp[0]; 
  // Amplitude(s) for diagram number 113
  FFV1_0(w[12], w[14], w[1], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[20] += -amp[0]; 
  VVVV1P0_1(w[0], w[1], w[4], cxtype(cIPC[4], cIPC[5]), 0., 0., w[12]); 
  VVVV3P0_1(w[0], w[1], w[4], cxtype(cIPC[4], cIPC[5]), 0., 0., w[24]); 
  VVVV4P0_1(w[0], w[1], w[4], cxtype(cIPC[4], cIPC[5]), 0., 0., w[21]); 
  // Amplitude(s) for diagram number 114
  VVV1_0(w[12], w[7], w[5], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[14] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[24], w[7], w[5], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[8] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[21], w[7], w[5], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += -cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  jamp[23] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 115
  FFV1_0(w[3], w[14], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += +amp[0]; 
  jamp[19] += -amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[23] += +amp[0]; 
  FFV1_0(w[3], w[14], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[19] += -amp[0]; 
  jamp[20] += +amp[0]; 
  jamp[21] += -amp[0]; 
  jamp[22] += +amp[0]; 
  FFV1_0(w[3], w[14], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[18] += -amp[0]; 
  jamp[20] += +amp[0]; 
  jamp[22] += +amp[0]; 
  jamp[23] += -amp[0]; 
  // Amplitude(s) for diagram number 116
  FFV1_0(w[13], w[2], w[12], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += +amp[0]; 
  jamp[2] += -amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[14] += +amp[0]; 
  FFV1_0(w[13], w[2], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[2] += -amp[0]; 
  jamp[6] += +amp[0]; 
  jamp[8] += -amp[0]; 
  jamp[12] += +amp[0]; 
  FFV1_0(w[13], w[2], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[0] += -amp[0]; 
  jamp[6] += +amp[0]; 
  jamp[12] += +amp[0]; 
  jamp[14] += -amp[0]; 
  VVVV1P0_1(w[0], w[1], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[21]); 
  VVVV3P0_1(w[0], w[1], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[13]); 
  VVVV4P0_1(w[0], w[1], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[24]); 
  // Amplitude(s) for diagram number 117
  VVV1_0(w[21], w[7], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[12] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[20] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[13], w[7], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[24], w[7], w[4], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += -cxtype(0, 1) * amp[0]; 
  jamp[7] += +cxtype(0, 1) * amp[0]; 
  jamp[12] += +cxtype(0, 1) * amp[0]; 
  jamp[14] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += +cxtype(0, 1) * amp[0]; 
  jamp[18] += +cxtype(0, 1) * amp[0]; 
  jamp[20] += -cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 118
  FFV1_0(w[3], w[11], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += +amp[0]; 
  jamp[13] += -amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[17] += +amp[0]; 
  FFV1_0(w[3], w[11], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[13] += -amp[0]; 
  jamp[14] += +amp[0]; 
  jamp[15] += -amp[0]; 
  jamp[16] += +amp[0]; 
  FFV1_0(w[3], w[11], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[12] += -amp[0]; 
  jamp[14] += +amp[0]; 
  jamp[16] += +amp[0]; 
  jamp[17] += -amp[0]; 
  // Amplitude(s) for diagram number 119
  FFV1_0(w[15], w[2], w[21], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += +amp[0]; 
  jamp[4] += -amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[20] += +amp[0]; 
  FFV1_0(w[15], w[2], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[4] += -amp[0]; 
  jamp[7] += +amp[0]; 
  jamp[10] += -amp[0]; 
  jamp[18] += +amp[0]; 
  FFV1_0(w[15], w[2], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[1] += -amp[0]; 
  jamp[7] += +amp[0]; 
  jamp[18] += +amp[0]; 
  jamp[20] += -amp[0]; 
  VVVV1P0_1(w[0], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[24]); 
  VVVV3P0_1(w[0], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[15]); 
  VVVV4P0_1(w[0], w[4], w[5], cxtype(cIPC[4], cIPC[5]), 0., 0., w[13]); 
  // Amplitude(s) for diagram number 120
  FFV1_0(w[3], w[9], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += +amp[0]; 
  jamp[7] += -amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[11] += +amp[0]; 
  FFV1_0(w[3], w[9], w[15], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[7] += -amp[0]; 
  jamp[8] += +amp[0]; 
  jamp[9] += -amp[0]; 
  jamp[10] += +amp[0]; 
  FFV1_0(w[3], w[9], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[6] += -amp[0]; 
  jamp[8] += +amp[0]; 
  jamp[10] += +amp[0]; 
  jamp[11] += -amp[0]; 
  // Amplitude(s) for diagram number 121
  FFV1_0(w[16], w[2], w[24], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += +amp[0]; 
  jamp[5] += -amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[22] += +amp[0]; 
  FFV1_0(w[16], w[2], w[15], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[5] += -amp[0]; 
  jamp[13] += +amp[0]; 
  jamp[16] += -amp[0]; 
  jamp[19] += +amp[0]; 
  FFV1_0(w[16], w[2], w[13], cxtype(cIPC[2], cIPC[3]), &amp[0]); 
  jamp[3] += -amp[0]; 
  jamp[13] += +amp[0]; 
  jamp[19] += +amp[0]; 
  jamp[22] += -amp[0]; 
  // Amplitude(s) for diagram number 122
  VVV1_0(w[24], w[1], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += +cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[22] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[15], w[1], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[7] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[16] += +cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  VVV1_0(w[13], w[1], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[6] += -cxtype(0, 1) * amp[0]; 
  jamp[8] += +cxtype(0, 1) * amp[0]; 
  jamp[10] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[13] += -cxtype(0, 1) * amp[0]; 
  jamp[19] += -cxtype(0, 1) * amp[0]; 
  jamp[22] += +cxtype(0, 1) * amp[0]; 
  // Amplitude(s) for diagram number 123
  VVV1_0(w[0], w[17], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += -cxtype(0, 1) * amp[0]; 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[5] += -cxtype(0, 1) * amp[0]; 
  jamp[9] += +cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[23] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[0], w[19], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[1] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[3] += +cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[11] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[17] += -cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  VVV1_0(w[0], w[8], w[7], cxtype(cIPC[0], cIPC[1]), &amp[0]); 
  jamp[0] += +cxtype(0, 1) * amp[0]; 
  jamp[2] += -cxtype(0, 1) * amp[0]; 
  jamp[4] += -cxtype(0, 1) * amp[0]; 
  jamp[5] += +cxtype(0, 1) * amp[0]; 
  jamp[9] += -cxtype(0, 1) * amp[0]; 
  jamp[15] += +cxtype(0, 1) * amp[0]; 
  jamp[21] += +cxtype(0, 1) * amp[0]; 
  jamp[23] += -cxtype(0, 1) * amp[0]; 
  // double CPPProcess::matrix_1_gg_ttxgg() {



  // Local variables

  // The color matrix;
  static const double denom[ncolor] = {54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54};
  static const double cf[ncolor][ncolor] = {{512, -64, -64, 8, 8, 80, -64, 8,
      8, -1, -1, -10, 8, -1, 80, -10, 71, 62, -1, -10, -10, 62, 62, -28}, {-64,
      512, 8, 80, -64, 8, 8, -64, -1, -10, 8, -1, -1, -10, -10, 62, 62, -28, 8,
      -1, 80, -10, 71, 62}, {-64, 8, 512, -64, 80, 8, 8, -1, 80, -10, 71, 62,
      -64, 8, 8, -1, -1, -10, -10, -1, 62, -28, -10, 62}, {8, 80, -64, 512, 8,
      -64, -1, -10, -10, 62, 62, -28, 8, -64, -1, -10, 8, -1, -1, 8, 71, 62,
      80, -10}, {8, -64, 80, 8, 512, -64, -1, 8, 71, 62, 80, -10, -10, -1, 62,
      -28, -10, 62, -64, 8, 8, -1, -1, -10}, {80, 8, 8, -64, -64, 512, -10, -1,
      62, -28, -10, 62, -1, 8, 71, 62, 80, -10, 8, -64, -1, -10, 8, -1}, {-64,
      8, 8, -1, -1, -10, 512, -64, -64, 8, 8, 80, 80, -10, 8, -1, 62, 71, -10,
      62, -1, -10, -28, 62}, {8, -64, -1, -10, 8, -1, -64, 512, 8, 80, -64, 8,
      -10, 62, -1, -10, -28, 62, 80, -10, 8, -1, 62, 71}, {8, -1, 80, -10, 71,
      62, -64, 8, 512, -64, 80, 8, 8, -1, -64, 8, -10, -1, 62, -28, -10, -1,
      62, -10}, {-1, -10, -10, 62, 62, -28, 8, 80, -64, 512, 8, -64, -1, -10,
      8, -64, -1, 8, 71, 62, -1, 8, -10, 80}, {-1, 8, 71, 62, 80, -10, 8, -64,
      80, 8, 512, -64, 62, -28, -10, -1, 62, -10, 8, -1, -64, 8, -10, -1},
      {-10, -1, 62, -28, -10, 62, 80, 8, 8, -64, -64, 512, 71, 62, -1, 8, -10,
      80, -1, -10, 8, -64, -1, 8}, {8, -1, -64, 8, -10, -1, 80, -10, 8, -1, 62,
      71, 512, -64, -64, 8, 8, 80, 62, -10, -28, 62, -1, -10}, {-1, -10, 8,
      -64, -1, 8, -10, 62, -1, -10, -28, 62, -64, 512, 8, 80, -64, 8, -10, 80,
      62, 71, 8, -1}, {80, -10, 8, -1, 62, 71, 8, -1, -64, 8, -10, -1, -64, 8,
      512, -64, 80, 8, -28, 62, 62, -10, -10, -1}, {-10, 62, -1, -10, -28, 62,
      -1, -10, 8, -64, -1, 8, 8, 80, -64, 512, 8, -64, 62, 71, -10, 80, -1, 8},
      {71, 62, -1, 8, -10, 80, 62, -28, -10, -1, 62, -10, 8, -64, 80, 8, 512,
      -64, -1, 8, -10, -1, -64, 8}, {62, -28, -10, -1, 62, -10, 71, 62, -1, 8,
      -10, 80, 80, 8, 8, -64, -64, 512, -10, -1, -1, 8, 8, -64}, {-1, 8, -10,
      -1, -64, 8, -10, 80, 62, 71, 8, -1, 62, -10, -28, 62, -1, -10, 512, -64,
      -64, 8, 8, 80}, {-10, -1, -1, 8, 8, -64, 62, -10, -28, 62, -1, -10, -10,
      80, 62, 71, 8, -1, -64, 512, 8, 80, -64, 8}, {-10, 80, 62, 71, 8, -1, -1,
      8, -10, -1, -64, 8, -28, 62, 62, -10, -10, -1, -64, 8, 512, -64, 80, 8},
      {62, -10, -28, 62, -1, -10, -10, -1, -1, 8, 8, -64, 62, 71, -10, 80, -1,
      8, 8, 80, -64, 512, 8, -64}, {62, 71, -10, 80, -1, 8, -28, 62, 62, -10,
      -10, -1, -1, 8, -10, -1, -64, 8, 8, -64, 80, 8, 512, -64}, {-28, 62, 62,
      -10, -10, -1, 62, 71, -10, 80, -1, 8, -10, -1, -1, 8, 8, -64, 80, 8, 8,
      -64, -64, 512}};


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
  static const int tHel[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1}, {-1, -1,
      -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1}, {-1, -1,
      -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1}, {-1, -1,
      -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1}, {-1, -1,
      1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1}, {-1, -1, 1,
      1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1, 1, -1, -1,
      -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1, 1, -1, -1,
      1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, 1,
      -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1, -1, -1, 1},
      {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1, -1, -1}, {-1,
      1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1}, {1, -1, -1,
      -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1}, {1, -1, -1,
      -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1, -1, -1, 1,
      1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, -1,
      1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1, -1, -1},
      {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1}, {1, 1,
      -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1, 1, -1,
      -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1, 1, 1,
      -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1, 1},
      {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1, 1,
      1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
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
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->mdl_MT); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 

  static cxtype tIPC[3] = {pars->GC_10, pars->GC_11, pars->GC_12}; 
  static double tIPD[2] = {pars->mdl_MT, pars->mdl_WT}; 

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
    const int denominators[1] = {512}; 

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


