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
// Process: g g > t t~ g g WEIGHTED<=4 @1

__constant__ int cHel[64][6]; 
// __constant__ double cmME[6]; value hardcoded now
// extern __constant__ int cPerm[4];
// 
__constant__ double cIPC[6];  // coupling ?
__constant__ double cIPD[2]; 


// Evaluate |M|^2 for each subprocess

__device__ void calculate_wavefunctions(int ihel, double local_mom[6][3],
    double &matrix)
{
  thrust::complex<double> amp[159]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[26][6]; 
  vxxxxx(local_mom[0], 0., cHel[ihel][0], -1, w[0]); 
  vxxxxx(local_mom[1], 0., cHel[ihel][1], -1, w[1]); 
  oxxxxx(local_mom[2], cIPD[0], cHel[ihel][2], +1, w[2]); 
  ixxxxx(local_mom[3], cIPD[0], cHel[ihel][3], -1, w[3]); 
  vxxxxx(local_mom[4], 0., cHel[ihel][4], +1, w[4]); 
  vxxxxx(local_mom[5], 0., cHel[ihel][5], +1, w[5]); 
  VVV1P0_1(w[0], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[6]);
  FFV1P0_3(w[3], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[7]);
  // Amplitude(s) for diagram number 1
  VVVV1_0(w[6], w[7], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[0]);
  VVVV3_0(w[6], w[7], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[1]);
  VVVV4_0(w[6], w[7], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[2]);
  VVV1P0_1(w[6], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[8]);
  // Amplitude(s) for diagram number 2
  VVV1_0(w[7], w[5], w[8], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[3]);
  VVV1P0_1(w[6], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[9]);
  // Amplitude(s) for diagram number 3
  VVV1_0(w[7], w[4], w[9], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[4]);
  VVV1P0_1(w[4], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[10]);
  // Amplitude(s) for diagram number 4
  VVV1_0(w[6], w[7], w[10], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[5]);
  FFV1_1(w[2], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[11]);
  FFV1_2(w[3], w[6], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  // Amplitude(s) for diagram number 5
  FFV1_0(w[12], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[6]);
  // Amplitude(s) for diagram number 6
  FFV1_0(w[3], w[11], w[9], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[7]);
  FFV1_2(w[3], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[13]);
  // Amplitude(s) for diagram number 7
  FFV1_0(w[13], w[11], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[8]);
  FFV1_1(w[2], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[14]);
  // Amplitude(s) for diagram number 8
  FFV1_0(w[12], w[14], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[9]);
  // Amplitude(s) for diagram number 9
  FFV1_0(w[3], w[14], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[10]);
  FFV1_2(w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[15]);
  // Amplitude(s) for diagram number 10
  FFV1_0(w[15], w[14], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[11]);
  FFV1_1(w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[16]);
  // Amplitude(s) for diagram number 11
  FFV1_0(w[15], w[16], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[12]);
  // Amplitude(s) for diagram number 12
  FFV1_0(w[15], w[2], w[9], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[13]);
  // Amplitude(s) for diagram number 13
  FFV1_0(w[13], w[16], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[14]);
  // Amplitude(s) for diagram number 14
  FFV1_0(w[13], w[2], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[15]);
  // Amplitude(s) for diagram number 15
  FFV1_0(w[3], w[16], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[16]);
  // Amplitude(s) for diagram number 16
  FFV1_0(w[12], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[17]);
  FFV1_1(w[2], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  FFV1_2(w[3], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[16]);
  FFV1_1(w[12], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[8]);
  // Amplitude(s) for diagram number 17
  FFV1_0(w[16], w[8], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[18]);
  FFV1_1(w[12], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[9]);
  // Amplitude(s) for diagram number 18
  FFV1_0(w[16], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[19]);
  // Amplitude(s) for diagram number 19
  FFV1_0(w[16], w[12], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[20]);
  VVV1P0_1(w[1], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[6]);
  FFV1P0_3(w[3], w[12], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[17]);
  // Amplitude(s) for diagram number 20
  VVV1_0(w[6], w[5], w[17], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[21]);
  // Amplitude(s) for diagram number 21
  FFV1_0(w[3], w[9], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[22]);
  // Amplitude(s) for diagram number 22
  FFV1_0(w[13], w[12], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[23]);
  VVV1P0_1(w[1], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[18]);
  // Amplitude(s) for diagram number 23
  VVV1_0(w[18], w[4], w[17], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[24]);
  // Amplitude(s) for diagram number 24
  FFV1_0(w[3], w[8], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[25]);
  // Amplitude(s) for diagram number 25
  FFV1_0(w[15], w[12], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[26]);
  FFV1_1(w[12], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[19]);
  // Amplitude(s) for diagram number 26
  FFV1_0(w[15], w[19], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[27]);
  // Amplitude(s) for diagram number 27
  FFV1_0(w[15], w[9], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[28]);
  // Amplitude(s) for diagram number 28
  FFV1_0(w[13], w[19], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[29]);
  // Amplitude(s) for diagram number 29
  FFV1_0(w[13], w[8], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[30]);
  // Amplitude(s) for diagram number 30
  FFV1_0(w[3], w[19], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[31]);
  // Amplitude(s) for diagram number 31
  VVV1_0(w[1], w[10], w[17], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[32]);
  VVVV1P0_1(w[1], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[17]);
  VVVV3P0_1(w[1], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[19]);
  VVVV4P0_1(w[1], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[8]);
  // Amplitude(s) for diagram number 32
  FFV1_0(w[3], w[12], w[17], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[33]);
  FFV1_0(w[3], w[12], w[19], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[34]);
  FFV1_0(w[3], w[12], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[35]);
  FFV1_2(w[3], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  FFV1_1(w[2], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[9]);
  FFV1_2(w[12], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[20]);
  // Amplitude(s) for diagram number 33
  FFV1_0(w[20], w[9], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[36]);
  FFV1_2(w[12], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[21]);
  // Amplitude(s) for diagram number 34
  FFV1_0(w[21], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[37]);
  // Amplitude(s) for diagram number 35
  FFV1_0(w[12], w[9], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[38]);
  FFV1P0_3(w[12], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 36
  VVV1_0(w[6], w[5], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[39]);
  // Amplitude(s) for diagram number 37
  FFV1_0(w[21], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[40]);
  // Amplitude(s) for diagram number 38
  FFV1_0(w[12], w[14], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[41]);
  // Amplitude(s) for diagram number 39
  VVV1_0(w[18], w[4], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[42]);
  // Amplitude(s) for diagram number 40
  FFV1_0(w[20], w[2], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[43]);
  // Amplitude(s) for diagram number 41
  FFV1_0(w[12], w[11], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[44]);
  FFV1_2(w[12], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[23]);
  // Amplitude(s) for diagram number 42
  FFV1_0(w[23], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[45]);
  // Amplitude(s) for diagram number 43
  FFV1_0(w[21], w[11], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[46]);
  // Amplitude(s) for diagram number 44
  FFV1_0(w[23], w[14], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[47]);
  // Amplitude(s) for diagram number 45
  FFV1_0(w[20], w[14], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[48]);
  // Amplitude(s) for diagram number 46
  FFV1_0(w[23], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[49]);
  // Amplitude(s) for diagram number 47
  VVV1_0(w[1], w[10], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[50]);
  // Amplitude(s) for diagram number 48
  FFV1_0(w[12], w[2], w[17], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[51]);
  FFV1_0(w[12], w[2], w[19], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[52]);
  FFV1_0(w[12], w[2], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[53]);
  VVV1P0_1(w[0], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[12]);
  FFV1_2(w[3], w[12], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[22]);
  // Amplitude(s) for diagram number 49
  FFV1_0(w[22], w[9], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[54]);
  VVV1P0_1(w[12], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[23]);
  // Amplitude(s) for diagram number 50
  FFV1_0(w[3], w[9], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[55]);
  // Amplitude(s) for diagram number 51
  FFV1_0(w[13], w[9], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[56]);
  FFV1_1(w[2], w[12], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[20]);
  // Amplitude(s) for diagram number 52
  FFV1_0(w[16], w[20], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[57]);
  // Amplitude(s) for diagram number 53
  FFV1_0(w[16], w[2], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[58]);
  // Amplitude(s) for diagram number 54
  FFV1_0(w[16], w[14], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[59]);
  // Amplitude(s) for diagram number 55
  FFV1_0(w[3], w[20], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[60]);
  // Amplitude(s) for diagram number 56
  FFV1_0(w[22], w[2], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[61]);
  // Amplitude(s) for diagram number 57
  VVV1_0(w[12], w[18], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[62]);
  // Amplitude(s) for diagram number 58
  VVVV1_0(w[12], w[1], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[63]);
  VVVV3_0(w[12], w[1], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[64]);
  VVVV4_0(w[12], w[1], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[65]);
  VVV1P0_1(w[12], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[21]);
  // Amplitude(s) for diagram number 59
  VVV1_0(w[7], w[5], w[21], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[66]);
  // Amplitude(s) for diagram number 60
  VVV1_0(w[1], w[7], w[23], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[67]);
  // Amplitude(s) for diagram number 61
  FFV1_0(w[3], w[14], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[68]);
  // Amplitude(s) for diagram number 62
  FFV1_0(w[22], w[14], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[69]);
  // Amplitude(s) for diagram number 63
  FFV1_0(w[13], w[2], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[70]);
  // Amplitude(s) for diagram number 64
  FFV1_0(w[13], w[20], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[71]);
  VVV1P0_1(w[0], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[20]);
  FFV1_2(w[3], w[20], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[21]);
  // Amplitude(s) for diagram number 65
  FFV1_0(w[21], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[72]);
  VVV1P0_1(w[20], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 66
  FFV1_0(w[3], w[9], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[73]);
  // Amplitude(s) for diagram number 67
  FFV1_0(w[15], w[9], w[20], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[74]);
  FFV1_1(w[2], w[20], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[23]);
  // Amplitude(s) for diagram number 68
  FFV1_0(w[16], w[23], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[75]);
  // Amplitude(s) for diagram number 69
  FFV1_0(w[16], w[2], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[76]);
  // Amplitude(s) for diagram number 70
  FFV1_0(w[16], w[11], w[20], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[77]);
  // Amplitude(s) for diagram number 71
  FFV1_0(w[3], w[23], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[78]);
  // Amplitude(s) for diagram number 72
  FFV1_0(w[21], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[79]);
  // Amplitude(s) for diagram number 73
  VVV1_0(w[20], w[6], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[80]);
  // Amplitude(s) for diagram number 74
  VVVV1_0(w[20], w[1], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[81]);
  VVVV3_0(w[20], w[1], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[82]);
  VVVV4_0(w[20], w[1], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[83]);
  VVV1P0_1(w[20], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[12]);
  // Amplitude(s) for diagram number 75
  VVV1_0(w[7], w[4], w[12], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[84]);
  // Amplitude(s) for diagram number 76
  VVV1_0(w[1], w[7], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[85]);
  // Amplitude(s) for diagram number 77
  FFV1_0(w[3], w[11], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[86]);
  // Amplitude(s) for diagram number 78
  FFV1_0(w[21], w[11], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[87]);
  // Amplitude(s) for diagram number 79
  FFV1_0(w[15], w[2], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[88]);
  // Amplitude(s) for diagram number 80
  FFV1_0(w[15], w[23], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[89]);
  FFV1_1(w[9], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[23]);
  // Amplitude(s) for diagram number 81
  FFV1_0(w[15], w[23], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[90]);
  FFV1_2(w[15], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  // Amplitude(s) for diagram number 82
  FFV1_0(w[12], w[9], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[91]);
  // Amplitude(s) for diagram number 83
  FFV1_0(w[13], w[23], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[92]);
  FFV1_2(w[13], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[21]);
  // Amplitude(s) for diagram number 84
  FFV1_0(w[21], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[93]);
  // Amplitude(s) for diagram number 85
  FFV1_0(w[3], w[23], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[94]);
  VVV1P0_1(w[0], w[10], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[23]);
  // Amplitude(s) for diagram number 86
  FFV1_0(w[3], w[9], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[95]);
  FFV1_2(w[16], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[22]);
  // Amplitude(s) for diagram number 87
  FFV1_0(w[22], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[96]);
  FFV1_1(w[11], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[20]);
  // Amplitude(s) for diagram number 88
  FFV1_0(w[16], w[20], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[97]);
  // Amplitude(s) for diagram number 89
  FFV1_0(w[22], w[14], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[98]);
  FFV1_1(w[14], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[24]);
  // Amplitude(s) for diagram number 90
  FFV1_0(w[16], w[24], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[99]);
  // Amplitude(s) for diagram number 91
  FFV1_0(w[22], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[100]);
  // Amplitude(s) for diagram number 92
  FFV1_0(w[16], w[2], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[101]);
  // Amplitude(s) for diagram number 93
  VVVV1_0(w[0], w[6], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[102]);
  VVVV3_0(w[0], w[6], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[103]);
  VVVV4_0(w[0], w[6], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[104]);
  VVV1P0_1(w[0], w[6], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 94
  VVV1_0(w[7], w[5], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[105]);
  VVV1P0_1(w[0], w[7], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[25]);
  // Amplitude(s) for diagram number 95
  VVV1_0(w[6], w[5], w[25], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[106]);
  // Amplitude(s) for diagram number 96
  FFV1_0(w[3], w[14], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[107]);
  // Amplitude(s) for diagram number 97
  FFV1_0(w[3], w[24], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[108]);
  // Amplitude(s) for diagram number 98
  FFV1_0(w[13], w[2], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[109]);
  // Amplitude(s) for diagram number 99
  FFV1_0(w[21], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[110]);
  // Amplitude(s) for diagram number 100
  VVVV1_0(w[0], w[18], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[111]);
  VVVV3_0(w[0], w[18], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[112]);
  VVVV4_0(w[0], w[18], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[113]);
  VVV1P0_1(w[0], w[18], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[6]);
  // Amplitude(s) for diagram number 101
  VVV1_0(w[7], w[4], w[6], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[114]);
  // Amplitude(s) for diagram number 102
  VVV1_0(w[18], w[4], w[25], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[115]);
  // Amplitude(s) for diagram number 103
  FFV1_0(w[3], w[11], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[116]);
  // Amplitude(s) for diagram number 104
  FFV1_0(w[3], w[20], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[117]);
  // Amplitude(s) for diagram number 105
  FFV1_0(w[15], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[118]);
  // Amplitude(s) for diagram number 106
  FFV1_0(w[12], w[2], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[119]);
  // Amplitude(s) for diagram number 107
  VVVV1_0(w[0], w[1], w[7], w[10], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[120]);
  VVVV3_0(w[0], w[1], w[7], w[10], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[121]);
  VVVV4_0(w[0], w[1], w[7], w[10], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[122]);
  // Amplitude(s) for diagram number 108
  VVV1_0(w[1], w[10], w[25], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[123]);
  // Amplitude(s) for diagram number 109
  VVV1_0(w[1], w[7], w[23], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[124]);
  // Amplitude(s) for diagram number 110
  FFV1_0(w[13], w[20], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[125]);
  // Amplitude(s) for diagram number 111
  FFV1_0(w[21], w[11], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[126]);
  // Amplitude(s) for diagram number 112
  FFV1_0(w[15], w[24], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[127]);
  // Amplitude(s) for diagram number 113
  FFV1_0(w[12], w[14], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[128]);
  VVVV1P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[12]);
  VVVV3P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[24]);
  VVVV4P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[21]);
  // Amplitude(s) for diagram number 114
  VVV1_0(w[12], w[7], w[5], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[129]);
  VVV1_0(w[24], w[7], w[5], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[130]);
  VVV1_0(w[21], w[7], w[5], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[131]);
  // Amplitude(s) for diagram number 115
  FFV1_0(w[3], w[14], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[132]);
  FFV1_0(w[3], w[14], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[133]);
  FFV1_0(w[3], w[14], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[134]);
  // Amplitude(s) for diagram number 116
  FFV1_0(w[13], w[2], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[135]);
  FFV1_0(w[13], w[2], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[136]);
  FFV1_0(w[13], w[2], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[137]);
  VVVV1P0_1(w[0], w[1], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[21]);
  VVVV3P0_1(w[0], w[1], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[13]);
  VVVV4P0_1(w[0], w[1], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[24]);
  // Amplitude(s) for diagram number 117
  VVV1_0(w[21], w[7], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[138]);
  VVV1_0(w[13], w[7], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[139]);
  VVV1_0(w[24], w[7], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[140]);
  // Amplitude(s) for diagram number 118
  FFV1_0(w[3], w[11], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[141]);
  FFV1_0(w[3], w[11], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[142]);
  FFV1_0(w[3], w[11], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[143]);
  // Amplitude(s) for diagram number 119
  FFV1_0(w[15], w[2], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[144]);
  FFV1_0(w[15], w[2], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[145]);
  FFV1_0(w[15], w[2], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[146]);
  VVVV1P0_1(w[0], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[24]);
  VVVV3P0_1(w[0], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[15]);
  VVVV4P0_1(w[0], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[13]);
  // Amplitude(s) for diagram number 120
  FFV1_0(w[3], w[9], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[147]);
  FFV1_0(w[3], w[9], w[15], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[148]);
  FFV1_0(w[3], w[9], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[149]);
  // Amplitude(s) for diagram number 121
  FFV1_0(w[16], w[2], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[150]);
  FFV1_0(w[16], w[2], w[15], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[151]);
  FFV1_0(w[16], w[2], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[152]);
  // Amplitude(s) for diagram number 122
  VVV1_0(w[24], w[1], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[153]);
  VVV1_0(w[15], w[1], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[154]);
  VVV1_0(w[13], w[1], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[155]);
  // Amplitude(s) for diagram number 123
  VVV1_0(w[0], w[17], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[156]);
  VVV1_0(w[0], w[19], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[157]);
  VVV1_0(w[0], w[8], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[158]);
  // double CPPProcess::matrix_1_gg_ttxgg() {
  int i, j; 
  // Local variables

  // const int ngraphs = 2;
  const int ncolor = 24; 
  thrust::complex<double> ztemp; 
  thrust::complex<double> jamp[ncolor]; 
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

  // Calculate color flows
  jamp[0] = +thrust::complex<double> (0, 1) * amp[0] + thrust::complex<double>
      (0, 1) * amp[1] + thrust::complex<double> (0, 1) * amp[3] +
      thrust::complex<double> (0, 1) * amp[5] + thrust::complex<double> (0, 1)
      * amp[14] + amp[15] + amp[16] + amp[21] + thrust::complex<double> (0, 1)
      * amp[23] - amp[29] + thrust::complex<double> (0, 1) * amp[31] + amp[32]
      + amp[33] - amp[35] + thrust::complex<double> (0, 1) * amp[102] +
      thrust::complex<double> (0, 1) * amp[103] + thrust::complex<double> (0,
      1) * amp[105] + thrust::complex<double> (0, 1) * amp[106] + amp[109] +
      thrust::complex<double> (0, 1) * amp[120] + thrust::complex<double> (0,
      1) * amp[121] + thrust::complex<double> (0, 1) * amp[123] +
      thrust::complex<double> (0, 1) * amp[129] - thrust::complex<double> (0,
      1) * amp[131] + amp[135] - amp[137] - thrust::complex<double> (0, 1) *
      amp[156] + thrust::complex<double> (0, 1) * amp[158];
  jamp[1] = -thrust::complex<double> (0, 1) * amp[0] + thrust::complex<double>
      (0, 1) * amp[2] + thrust::complex<double> (0, 1) * amp[4] -
      thrust::complex<double> (0, 1) * amp[5] + thrust::complex<double> (0, 1)
      * amp[12] + amp[13] - amp[16] + amp[24] + thrust::complex<double> (0, 1)
      * amp[26] - amp[27] - thrust::complex<double> (0, 1) * amp[31] - amp[32]
      - amp[33] - amp[34] + thrust::complex<double> (0, 1) * amp[111] +
      thrust::complex<double> (0, 1) * amp[112] + thrust::complex<double> (0,
      1) * amp[114] + thrust::complex<double> (0, 1) * amp[115] + amp[118] -
      thrust::complex<double> (0, 1) * amp[120] - thrust::complex<double> (0,
      1) * amp[121] - thrust::complex<double> (0, 1) * amp[123] +
      thrust::complex<double> (0, 1) * amp[138] - thrust::complex<double> (0,
      1) * amp[140] + amp[144] - amp[146] + thrust::complex<double> (0, 1) *
      amp[156] + thrust::complex<double> (0, 1) * amp[157];
  jamp[2] = -amp[21] - thrust::complex<double> (0, 1) * amp[23] - amp[24] +
      thrust::complex<double> (0, 1) * amp[25] - amp[30] + amp[34] + amp[35] +
      amp[60] - thrust::complex<double> (0, 1) * amp[62] +
      thrust::complex<double> (0, 1) * amp[63] + thrust::complex<double> (0, 1)
      * amp[64] + thrust::complex<double> (0, 1) * amp[66] + amp[70] +
      thrust::complex<double> (0, 1) * amp[71] - thrust::complex<double> (0, 1)
      * amp[102] - thrust::complex<double> (0, 1) * amp[103] -
      thrust::complex<double> (0, 1) * amp[105] - thrust::complex<double> (0,
      1) * amp[106] - amp[109] - thrust::complex<double> (0, 1) * amp[112] -
      thrust::complex<double> (0, 1) * amp[113] - thrust::complex<double> (0,
      1) * amp[115] - thrust::complex<double> (0, 1) * amp[129] -
      thrust::complex<double> (0, 1) * amp[130] - amp[135] - amp[136] -
      thrust::complex<double> (0, 1) * amp[157] - thrust::complex<double> (0,
      1) * amp[158];
  jamp[3] = -amp[18] + thrust::complex<double> (0, 1) * amp[20] + amp[24] -
      thrust::complex<double> (0, 1) * amp[25] - amp[32] - amp[33] - amp[34] +
      thrust::complex<double> (0, 1) * amp[57] + amp[58] - amp[60] +
      thrust::complex<double> (0, 1) * amp[62] - thrust::complex<double> (0, 1)
      * amp[64] - thrust::complex<double> (0, 1) * amp[65] -
      thrust::complex<double> (0, 1) * amp[67] + amp[101] +
      thrust::complex<double> (0, 1) * amp[112] + thrust::complex<double> (0,
      1) * amp[113] + thrust::complex<double> (0, 1) * amp[115] -
      thrust::complex<double> (0, 1) * amp[121] - thrust::complex<double> (0,
      1) * amp[122] - thrust::complex<double> (0, 1) * amp[123] -
      thrust::complex<double> (0, 1) * amp[124] + amp[150] - amp[152] -
      thrust::complex<double> (0, 1) * amp[153] + thrust::complex<double> (0,
      1) * amp[155] + thrust::complex<double> (0, 1) * amp[156] +
      thrust::complex<double> (0, 1) * amp[157];
  jamp[4] = -amp[21] + thrust::complex<double> (0, 1) * amp[22] - amp[24] -
      thrust::complex<double> (0, 1) * amp[26] - amp[28] + amp[34] + amp[35] +
      amp[78] - thrust::complex<double> (0, 1) * amp[80] +
      thrust::complex<double> (0, 1) * amp[81] + thrust::complex<double> (0, 1)
      * amp[82] + thrust::complex<double> (0, 1) * amp[84] + amp[88] +
      thrust::complex<double> (0, 1) * amp[89] - thrust::complex<double> (0, 1)
      * amp[103] - thrust::complex<double> (0, 1) * amp[104] -
      thrust::complex<double> (0, 1) * amp[106] - thrust::complex<double> (0,
      1) * amp[111] - thrust::complex<double> (0, 1) * amp[112] -
      thrust::complex<double> (0, 1) * amp[114] - thrust::complex<double> (0,
      1) * amp[115] - amp[118] - thrust::complex<double> (0, 1) * amp[138] -
      thrust::complex<double> (0, 1) * amp[139] - amp[144] - amp[145] -
      thrust::complex<double> (0, 1) * amp[157] - thrust::complex<double> (0,
      1) * amp[158];
  jamp[5] = -amp[19] - thrust::complex<double> (0, 1) * amp[20] + amp[21] -
      thrust::complex<double> (0, 1) * amp[22] + amp[32] + amp[33] - amp[35] +
      thrust::complex<double> (0, 1) * amp[75] + amp[76] - amp[78] +
      thrust::complex<double> (0, 1) * amp[80] - thrust::complex<double> (0, 1)
      * amp[82] - thrust::complex<double> (0, 1) * amp[83] -
      thrust::complex<double> (0, 1) * amp[85] - amp[101] +
      thrust::complex<double> (0, 1) * amp[103] + thrust::complex<double> (0,
      1) * amp[104] + thrust::complex<double> (0, 1) * amp[106] +
      thrust::complex<double> (0, 1) * amp[121] + thrust::complex<double> (0,
      1) * amp[122] + thrust::complex<double> (0, 1) * amp[123] +
      thrust::complex<double> (0, 1) * amp[124] - amp[150] - amp[151] +
      thrust::complex<double> (0, 1) * amp[153] + thrust::complex<double> (0,
      1) * amp[154] - thrust::complex<double> (0, 1) * amp[156] +
      thrust::complex<double> (0, 1) * amp[158];
  jamp[6] = -thrust::complex<double> (0, 1) * amp[0] - thrust::complex<double>
      (0, 1) * amp[1] - thrust::complex<double> (0, 1) * amp[3] -
      thrust::complex<double> (0, 1) * amp[5] - thrust::complex<double> (0, 1)
      * amp[14] - amp[15] - amp[16] + amp[55] + thrust::complex<double> (0, 1)
      * amp[56] - thrust::complex<double> (0, 1) * amp[63] +
      thrust::complex<double> (0, 1) * amp[65] - thrust::complex<double> (0, 1)
      * amp[66] + thrust::complex<double> (0, 1) * amp[67] - amp[70] - amp[92]
      + thrust::complex<double> (0, 1) * amp[94] + amp[95] -
      thrust::complex<double> (0, 1) * amp[120] + thrust::complex<double> (0,
      1) * amp[122] + thrust::complex<double> (0, 1) * amp[124] +
      thrust::complex<double> (0, 1) * amp[130] + thrust::complex<double> (0,
      1) * amp[131] + amp[136] + amp[137] + amp[147] - amp[149] +
      thrust::complex<double> (0, 1) * amp[153] - thrust::complex<double> (0,
      1) * amp[155];
  jamp[7] = +thrust::complex<double> (0, 1) * amp[0] - thrust::complex<double>
      (0, 1) * amp[2] - thrust::complex<double> (0, 1) * amp[4] +
      thrust::complex<double> (0, 1) * amp[5] - thrust::complex<double> (0, 1)
      * amp[12] - amp[13] + amp[16] + amp[73] + thrust::complex<double> (0, 1)
      * amp[74] - thrust::complex<double> (0, 1) * amp[81] +
      thrust::complex<double> (0, 1) * amp[83] - thrust::complex<double> (0, 1)
      * amp[84] + thrust::complex<double> (0, 1) * amp[85] - amp[88] - amp[90]
      - thrust::complex<double> (0, 1) * amp[94] - amp[95] +
      thrust::complex<double> (0, 1) * amp[120] - thrust::complex<double> (0,
      1) * amp[122] - thrust::complex<double> (0, 1) * amp[124] +
      thrust::complex<double> (0, 1) * amp[139] + thrust::complex<double> (0,
      1) * amp[140] + amp[145] + amp[146] - amp[147] - amp[148] -
      thrust::complex<double> (0, 1) * amp[153] - thrust::complex<double> (0,
      1) * amp[154];
  jamp[8] = -amp[55] - thrust::complex<double> (0, 1) * amp[56] +
      thrust::complex<double> (0, 1) * amp[63] - thrust::complex<double> (0, 1)
      * amp[65] + thrust::complex<double> (0, 1) * amp[66] -
      thrust::complex<double> (0, 1) * amp[67] + amp[70] +
      thrust::complex<double> (0, 1) * amp[72] - amp[73] + amp[79] +
      thrust::complex<double> (0, 1) * amp[80] - thrust::complex<double> (0, 1)
      * amp[82] - thrust::complex<double> (0, 1) * amp[83] -
      thrust::complex<double> (0, 1) * amp[85] - amp[93] -
      thrust::complex<double> (0, 1) * amp[102] + thrust::complex<double> (0,
      1) * amp[104] - thrust::complex<double> (0, 1) * amp[105] - amp[109] +
      thrust::complex<double> (0, 1) * amp[110] - thrust::complex<double> (0,
      1) * amp[129] - thrust::complex<double> (0, 1) * amp[130] - amp[135] -
      amp[136] + amp[148] + amp[149] + thrust::complex<double> (0, 1) *
      amp[154] + thrust::complex<double> (0, 1) * amp[155];
  jamp[9] = -amp[37] + thrust::complex<double> (0, 1) * amp[38] + amp[39] +
      thrust::complex<double> (0, 1) * amp[40] + amp[50] + amp[51] - amp[53] -
      thrust::complex<double> (0, 1) * amp[72] + amp[73] - amp[79] -
      thrust::complex<double> (0, 1) * amp[80] + thrust::complex<double> (0, 1)
      * amp[82] + thrust::complex<double> (0, 1) * amp[83] +
      thrust::complex<double> (0, 1) * amp[85] - amp[95] -
      thrust::complex<double> (0, 1) * amp[103] - thrust::complex<double> (0,
      1) * amp[104] - thrust::complex<double> (0, 1) * amp[106] -
      thrust::complex<double> (0, 1) * amp[121] - thrust::complex<double> (0,
      1) * amp[122] - thrust::complex<double> (0, 1) * amp[123] -
      thrust::complex<double> (0, 1) * amp[124] - amp[147] - amp[148] -
      thrust::complex<double> (0, 1) * amp[153] - thrust::complex<double> (0,
      1) * amp[154] + thrust::complex<double> (0, 1) * amp[156] -
      thrust::complex<double> (0, 1) * amp[158];
  jamp[10] = +thrust::complex<double> (0, 1) * amp[54] - amp[55] + amp[61] +
      thrust::complex<double> (0, 1) * amp[62] - thrust::complex<double> (0, 1)
      * amp[64] - thrust::complex<double> (0, 1) * amp[65] -
      thrust::complex<double> (0, 1) * amp[67] - amp[73] -
      thrust::complex<double> (0, 1) * amp[74] + thrust::complex<double> (0, 1)
      * amp[81] - thrust::complex<double> (0, 1) * amp[83] +
      thrust::complex<double> (0, 1) * amp[84] - thrust::complex<double> (0, 1)
      * amp[85] + amp[88] - amp[91] - thrust::complex<double> (0, 1) * amp[111]
      + thrust::complex<double> (0, 1) * amp[113] - thrust::complex<double> (0,
      1) * amp[114] - amp[118] + thrust::complex<double> (0, 1) * amp[119] -
      thrust::complex<double> (0, 1) * amp[138] - thrust::complex<double> (0,
      1) * amp[139] - amp[144] - amp[145] + amp[148] + amp[149] +
      thrust::complex<double> (0, 1) * amp[154] + thrust::complex<double> (0,
      1) * amp[155];
  jamp[11] = -amp[36] - thrust::complex<double> (0, 1) * amp[38] + amp[42] +
      thrust::complex<double> (0, 1) * amp[43] - amp[50] - amp[51] - amp[52] -
      thrust::complex<double> (0, 1) * amp[54] + amp[55] - amp[61] -
      thrust::complex<double> (0, 1) * amp[62] + thrust::complex<double> (0, 1)
      * amp[64] + thrust::complex<double> (0, 1) * amp[65] +
      thrust::complex<double> (0, 1) * amp[67] + amp[95] -
      thrust::complex<double> (0, 1) * amp[112] - thrust::complex<double> (0,
      1) * amp[113] - thrust::complex<double> (0, 1) * amp[115] +
      thrust::complex<double> (0, 1) * amp[121] + thrust::complex<double> (0,
      1) * amp[122] + thrust::complex<double> (0, 1) * amp[123] +
      thrust::complex<double> (0, 1) * amp[124] + amp[147] - amp[149] +
      thrust::complex<double> (0, 1) * amp[153] - thrust::complex<double> (0,
      1) * amp[155] - thrust::complex<double> (0, 1) * amp[156] -
      thrust::complex<double> (0, 1) * amp[157];
  jamp[12] = -thrust::complex<double> (0, 1) * amp[1] - thrust::complex<double>
      (0, 1) * amp[2] - thrust::complex<double> (0, 1) * amp[3] -
      thrust::complex<double> (0, 1) * amp[4] + amp[7] +
      thrust::complex<double> (0, 1) * amp[8] - amp[15] - amp[60] +
      thrust::complex<double> (0, 1) * amp[62] - thrust::complex<double> (0, 1)
      * amp[63] - thrust::complex<double> (0, 1) * amp[64] -
      thrust::complex<double> (0, 1) * amp[66] - amp[70] -
      thrust::complex<double> (0, 1) * amp[71] - thrust::complex<double> (0, 1)
      * amp[111] + thrust::complex<double> (0, 1) * amp[113] -
      thrust::complex<double> (0, 1) * amp[114] + amp[116] +
      thrust::complex<double> (0, 1) * amp[117] - amp[125] +
      thrust::complex<double> (0, 1) * amp[130] + thrust::complex<double> (0,
      1) * amp[131] + amp[136] + amp[137] - thrust::complex<double> (0, 1) *
      amp[138] + thrust::complex<double> (0, 1) * amp[140] + amp[141] -
      amp[143];
  jamp[13] = -thrust::complex<double> (0, 1) * amp[57] - amp[58] + amp[60] -
      thrust::complex<double> (0, 1) * amp[62] + thrust::complex<double> (0, 1)
      * amp[64] + thrust::complex<double> (0, 1) * amp[65] +
      thrust::complex<double> (0, 1) * amp[67] - amp[76] +
      thrust::complex<double> (0, 1) * amp[77] - thrust::complex<double> (0, 1)
      * amp[81] + thrust::complex<double> (0, 1) * amp[83] -
      thrust::complex<double> (0, 1) * amp[84] + thrust::complex<double> (0, 1)
      * amp[85] + amp[86] - amp[97] + thrust::complex<double> (0, 1) * amp[111]
      - thrust::complex<double> (0, 1) * amp[113] + thrust::complex<double> (0,
      1) * amp[114] - amp[116] - thrust::complex<double> (0, 1) * amp[117] +
      thrust::complex<double> (0, 1) * amp[138] + thrust::complex<double> (0,
      1) * amp[139] - amp[141] - amp[142] + amp[151] + amp[152] -
      thrust::complex<double> (0, 1) * amp[154] - thrust::complex<double> (0,
      1) * amp[155];
  jamp[14] = +thrust::complex<double> (0, 1) * amp[1] + thrust::complex<double>
      (0, 1) * amp[2] + thrust::complex<double> (0, 1) * amp[3] +
      thrust::complex<double> (0, 1) * amp[4] - amp[7] -
      thrust::complex<double> (0, 1) * amp[8] + amp[15] - amp[79] -
      thrust::complex<double> (0, 1) * amp[80] + thrust::complex<double> (0, 1)
      * amp[81] + thrust::complex<double> (0, 1) * amp[82] +
      thrust::complex<double> (0, 1) * amp[84] - amp[86] +
      thrust::complex<double> (0, 1) * amp[87] + thrust::complex<double> (0, 1)
      * amp[102] - thrust::complex<double> (0, 1) * amp[104] +
      thrust::complex<double> (0, 1) * amp[105] + amp[109] -
      thrust::complex<double> (0, 1) * amp[110] - amp[126] +
      thrust::complex<double> (0, 1) * amp[129] - thrust::complex<double> (0,
      1) * amp[131] + amp[135] - amp[137] - thrust::complex<double> (0, 1) *
      amp[139] - thrust::complex<double> (0, 1) * amp[140] + amp[142] +
      amp[143];
  jamp[15] = -amp[39] - thrust::complex<double> (0, 1) * amp[40] - amp[42] +
      thrust::complex<double> (0, 1) * amp[44] - amp[46] + amp[52] + amp[53] +
      amp[79] + thrust::complex<double> (0, 1) * amp[80] -
      thrust::complex<double> (0, 1) * amp[81] - thrust::complex<double> (0, 1)
      * amp[82] - thrust::complex<double> (0, 1) * amp[84] + amp[86] -
      thrust::complex<double> (0, 1) * amp[87] + thrust::complex<double> (0, 1)
      * amp[103] + thrust::complex<double> (0, 1) * amp[104] +
      thrust::complex<double> (0, 1) * amp[106] + thrust::complex<double> (0,
      1) * amp[111] + thrust::complex<double> (0, 1) * amp[112] +
      thrust::complex<double> (0, 1) * amp[114] + thrust::complex<double> (0,
      1) * amp[115] - amp[116] + thrust::complex<double> (0, 1) * amp[138] +
      thrust::complex<double> (0, 1) * amp[139] - amp[141] - amp[142] +
      thrust::complex<double> (0, 1) * amp[157] + thrust::complex<double> (0,
      1) * amp[158];
  jamp[16] = -thrust::complex<double> (0, 1) * amp[0] + thrust::complex<double>
      (0, 1) * amp[2] + thrust::complex<double> (0, 1) * amp[4] -
      thrust::complex<double> (0, 1) * amp[5] + thrust::complex<double> (0, 1)
      * amp[6] - amp[7] + amp[17] + amp[76] - thrust::complex<double> (0, 1) *
      amp[77] + thrust::complex<double> (0, 1) * amp[81] -
      thrust::complex<double> (0, 1) * amp[83] + thrust::complex<double> (0, 1)
      * amp[84] - thrust::complex<double> (0, 1) * amp[85] - amp[86] - amp[96]
      + thrust::complex<double> (0, 1) * amp[100] - amp[101] -
      thrust::complex<double> (0, 1) * amp[120] + thrust::complex<double> (0,
      1) * amp[122] + thrust::complex<double> (0, 1) * amp[124] -
      thrust::complex<double> (0, 1) * amp[139] - thrust::complex<double> (0,
      1) * amp[140] + amp[142] + amp[143] - amp[150] - amp[151] +
      thrust::complex<double> (0, 1) * amp[153] + thrust::complex<double> (0,
      1) * amp[154];
  jamp[17] = +thrust::complex<double> (0, 1) * amp[0] - thrust::complex<double>
      (0, 1) * amp[2] - thrust::complex<double> (0, 1) * amp[4] +
      thrust::complex<double> (0, 1) * amp[5] - thrust::complex<double> (0, 1)
      * amp[6] + amp[7] - amp[17] + amp[42] - thrust::complex<double> (0, 1) *
      amp[44] - amp[45] + thrust::complex<double> (0, 1) * amp[49] - amp[50] -
      amp[51] - amp[52] - thrust::complex<double> (0, 1) * amp[111] -
      thrust::complex<double> (0, 1) * amp[112] - thrust::complex<double> (0,
      1) * amp[114] - thrust::complex<double> (0, 1) * amp[115] + amp[116] +
      thrust::complex<double> (0, 1) * amp[120] + thrust::complex<double> (0,
      1) * amp[121] + thrust::complex<double> (0, 1) * amp[123] -
      thrust::complex<double> (0, 1) * amp[138] + thrust::complex<double> (0,
      1) * amp[140] + amp[141] - amp[143] - thrust::complex<double> (0, 1) *
      amp[156] - thrust::complex<double> (0, 1) * amp[157];
  jamp[18] = -thrust::complex<double> (0, 1) * amp[1] - thrust::complex<double>
      (0, 1) * amp[2] - thrust::complex<double> (0, 1) * amp[3] -
      thrust::complex<double> (0, 1) * amp[4] + amp[10] +
      thrust::complex<double> (0, 1) * amp[11] - amp[13] - amp[78] +
      thrust::complex<double> (0, 1) * amp[80] - thrust::complex<double> (0, 1)
      * amp[81] - thrust::complex<double> (0, 1) * amp[82] -
      thrust::complex<double> (0, 1) * amp[84] - amp[88] -
      thrust::complex<double> (0, 1) * amp[89] - thrust::complex<double> (0, 1)
      * amp[102] + thrust::complex<double> (0, 1) * amp[104] -
      thrust::complex<double> (0, 1) * amp[105] + amp[107] +
      thrust::complex<double> (0, 1) * amp[108] - amp[127] -
      thrust::complex<double> (0, 1) * amp[129] + thrust::complex<double> (0,
      1) * amp[131] + amp[132] - amp[134] + thrust::complex<double> (0, 1) *
      amp[139] + thrust::complex<double> (0, 1) * amp[140] + amp[145] +
      amp[146];
  jamp[19] = -amp[58] + thrust::complex<double> (0, 1) * amp[59] -
      thrust::complex<double> (0, 1) * amp[63] + thrust::complex<double> (0, 1)
      * amp[65] - thrust::complex<double> (0, 1) * amp[66] +
      thrust::complex<double> (0, 1) * amp[67] + amp[68] -
      thrust::complex<double> (0, 1) * amp[75] - amp[76] + amp[78] -
      thrust::complex<double> (0, 1) * amp[80] + thrust::complex<double> (0, 1)
      * amp[82] + thrust::complex<double> (0, 1) * amp[83] +
      thrust::complex<double> (0, 1) * amp[85] - amp[99] +
      thrust::complex<double> (0, 1) * amp[102] - thrust::complex<double> (0,
      1) * amp[104] + thrust::complex<double> (0, 1) * amp[105] - amp[107] -
      thrust::complex<double> (0, 1) * amp[108] + thrust::complex<double> (0,
      1) * amp[129] + thrust::complex<double> (0, 1) * amp[130] - amp[132] -
      amp[133] + amp[151] + amp[152] - thrust::complex<double> (0, 1) *
      amp[154] - thrust::complex<double> (0, 1) * amp[155];
  jamp[20] = +thrust::complex<double> (0, 1) * amp[1] + thrust::complex<double>
      (0, 1) * amp[2] + thrust::complex<double> (0, 1) * amp[3] +
      thrust::complex<double> (0, 1) * amp[4] - amp[10] -
      thrust::complex<double> (0, 1) * amp[11] + amp[13] - amp[61] -
      thrust::complex<double> (0, 1) * amp[62] + thrust::complex<double> (0, 1)
      * amp[63] + thrust::complex<double> (0, 1) * amp[64] +
      thrust::complex<double> (0, 1) * amp[66] - amp[68] +
      thrust::complex<double> (0, 1) * amp[69] + thrust::complex<double> (0, 1)
      * amp[111] - thrust::complex<double> (0, 1) * amp[113] +
      thrust::complex<double> (0, 1) * amp[114] + amp[118] -
      thrust::complex<double> (0, 1) * amp[119] - amp[128] -
      thrust::complex<double> (0, 1) * amp[130] - thrust::complex<double> (0,
      1) * amp[131] + amp[133] + amp[134] + thrust::complex<double> (0, 1) *
      amp[138] - thrust::complex<double> (0, 1) * amp[140] + amp[144] -
      amp[146];
  jamp[21] = -amp[39] + thrust::complex<double> (0, 1) * amp[41] - amp[42] -
      thrust::complex<double> (0, 1) * amp[43] - amp[48] + amp[52] + amp[53] +
      amp[61] + thrust::complex<double> (0, 1) * amp[62] -
      thrust::complex<double> (0, 1) * amp[63] - thrust::complex<double> (0, 1)
      * amp[64] - thrust::complex<double> (0, 1) * amp[66] + amp[68] -
      thrust::complex<double> (0, 1) * amp[69] + thrust::complex<double> (0, 1)
      * amp[102] + thrust::complex<double> (0, 1) * amp[103] +
      thrust::complex<double> (0, 1) * amp[105] + thrust::complex<double> (0,
      1) * amp[106] - amp[107] + thrust::complex<double> (0, 1) * amp[112] +
      thrust::complex<double> (0, 1) * amp[113] + thrust::complex<double> (0,
      1) * amp[115] + thrust::complex<double> (0, 1) * amp[129] +
      thrust::complex<double> (0, 1) * amp[130] - amp[132] - amp[133] +
      thrust::complex<double> (0, 1) * amp[157] + thrust::complex<double> (0,
      1) * amp[158];
  jamp[22] = +thrust::complex<double> (0, 1) * amp[0] + thrust::complex<double>
      (0, 1) * amp[1] + thrust::complex<double> (0, 1) * amp[3] +
      thrust::complex<double> (0, 1) * amp[5] + thrust::complex<double> (0, 1)
      * amp[9] - amp[10] - amp[17] + amp[58] - thrust::complex<double> (0, 1) *
      amp[59] + thrust::complex<double> (0, 1) * amp[63] -
      thrust::complex<double> (0, 1) * amp[65] + thrust::complex<double> (0, 1)
      * amp[66] - thrust::complex<double> (0, 1) * amp[67] - amp[68] - amp[98]
      - thrust::complex<double> (0, 1) * amp[100] + amp[101] +
      thrust::complex<double> (0, 1) * amp[120] - thrust::complex<double> (0,
      1) * amp[122] - thrust::complex<double> (0, 1) * amp[124] -
      thrust::complex<double> (0, 1) * amp[130] - thrust::complex<double> (0,
      1) * amp[131] + amp[133] + amp[134] + amp[150] - amp[152] -
      thrust::complex<double> (0, 1) * amp[153] + thrust::complex<double> (0,
      1) * amp[155];
  jamp[23] = -thrust::complex<double> (0, 1) * amp[0] - thrust::complex<double>
      (0, 1) * amp[1] - thrust::complex<double> (0, 1) * amp[3] -
      thrust::complex<double> (0, 1) * amp[5] - thrust::complex<double> (0, 1)
      * amp[9] + amp[10] + amp[17] + amp[39] - thrust::complex<double> (0, 1) *
      amp[41] - amp[47] - thrust::complex<double> (0, 1) * amp[49] + amp[50] +
      amp[51] - amp[53] - thrust::complex<double> (0, 1) * amp[102] -
      thrust::complex<double> (0, 1) * amp[103] - thrust::complex<double> (0,
      1) * amp[105] - thrust::complex<double> (0, 1) * amp[106] + amp[107] -
      thrust::complex<double> (0, 1) * amp[120] - thrust::complex<double> (0,
      1) * amp[121] - thrust::complex<double> (0, 1) * amp[123] -
      thrust::complex<double> (0, 1) * amp[129] + thrust::complex<double> (0,
      1) * amp[131] + amp[132] - amp[134] + thrust::complex<double> (0, 1) *
      amp[156] - thrust::complex<double> (0, 1) * amp[158];

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
  mME.push_back(pars->ZERO); 
  static thrust::complex<double> tIPC[3] = {pars->GC_10, pars->GC_11,
      pars->GC_12};
  static double tIPD[2] = {pars->mdl_MT, pars->mdl_WT}; 

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
  // size_t slicePitch = dpt * 6;

  // char *dps = devPtr + dim * slicePitch;
  double matrix_element[nprocesses]; 

  thrust::complex<double> amp[159]; 

  double local_m[6][3]; 
  int DIM = blockDim.x * gridDim.x; 
  // for (int i=0; i<20;i++){
  // printf(" %f ", allmomenta[i]);
  // }
  // printf("\n");
  // printf("DIM is %i/%i\n", tid, DIM);
  for (int i = 0; i < 6; i++ )
  {
    for (int j = 0; j < 3; j++ )
    {
      local_m[i][j] = allmomenta[i * 3 * DIM + j * DIM + tid]; 
      // printf(" %f ", local_m[i][j]);
    }
    // printf("\n");
  }


  // Local variables and constants
  const int ncomb = 64; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  // double t[1];
  // Helicities for the process
  // static const int helicities[ncomb][nexternal] =
  // {{-1,-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1,1},{-1,-1,-1,-1,1,-1},{-1,-1,-1,-1,1,1
  // },{-1,-1,-1,1,-1,-1},{-1,-1,-1,1,-1,1},{-1,-1,-1,1,1,-1},{-1,-1,-1,1,1,1},{
  // -1,-1,1,-1,-1,-1},{-1,-1,1,-1,-1,1},{-1,-1,1,-1,1,-1},{-1,-1,1,-1,1,1},{-1,
  // -1,1,1,-1,-1},{-1,-1,1,1,-1,1},{-1,-1,1,1,1,-1},{-1,-1,1,1,1,1},{-1,1,-1,-1
  // ,-1,-1},{-1,1,-1,-1,-1,1},{-1,1,-1,-1,1,-1},{-1,1,-1,-1,1,1},{-1,1,-1,1,-1,
  // -1},{-1,1,-1,1,-1,1},{-1,1,-1,1,1,-1},{-1,1,-1,1,1,1},{-1,1,1,-1,-1,-1},{-1
  // ,1,1,-1,-1,1},{-1,1,1,-1,1,-1},{-1,1,1,-1,1,1},{-1,1,1,1,-1,-1},{-1,1,1,1,-
  // 1,1},{-1,1,1,1,1,-1},{-1,1,1,1,1,1},{1,-1,-1,-1,-1,-1},{1,-1,-1,-1,-1,1},{1
  // ,-1,-1,-1,1,-1},{1,-1,-1,-1,1,1},{1,-1,-1,1,-1,-1},{1,-1,-1,1,-1,1},{1,-1,-
  // 1,1,1,-1},{1,-1,-1,1,1,1},{1,-1,1,-1,-1,-1},{1,-1,1,-1,-1,1},{1,-1,1,-1,1,-
  // 1},{1,-1,1,-1,1,1},{1,-1,1,1,-1,-1},{1,-1,1,1,-1,1},{1,-1,1,1,1,-1},{1,-1,1
  // ,1,1,1},{1,1,-1,-1,-1,-1},{1,1,-1,-1,-1,1},{1,1,-1,-1,1,-1},{1,1,-1,-1,1,1}
  // ,{1,1,-1,1,-1,-1},{1,1,-1,1,-1,1},{1,1,-1,1,1,-1},{1,1,-1,1,1,1},{1,1,1,-1,
  // -1,-1},{1,1,1,-1,-1,1},{1,1,1,-1,1,-1},{1,1,1,-1,1,1},{1,1,1,1,-1,-1},{1,1,
  // 1,1,-1,1},{1,1,1,1,1,-1},{1,1,1,1,1,1}};
  // Denominators: spins, colors and identical particles
  const int denominators[1] = {512}; 


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


