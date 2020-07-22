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

__device__ void calculate_wavefunctions(int ihel, char * dps, size_t dpt,
    double &matrix)
{
  thrust::complex<double> amp[159]; 
  // Calculate wavefunctions for all processes
  thrust::complex<double> w[26][6]; 
  vxxxxx((double * )(dps + 0 * dpt), 0., cHel[ihel][0], -1, w[0]); 
  vxxxxx((double * )(dps + 1 * dpt), 0., cHel[ihel][1], -1, w[1]); 
  oxxxxx((double * )(dps + 2 * dpt), cIPD[0], cHel[ihel][2], +1, w[2]); 
  ixxxxx((double * )(dps + 3 * dpt), cIPD[0], cHel[ihel][3], -1, w[3]); 
  vxxxxx((double * )(dps + 4 * dpt), 0., cHel[ihel][4], +1, w[4]); 
  vxxxxx((double * )(dps + 5 * dpt), 0., cHel[ihel][5], +1, w[5]); 
  VVV1P0_1(w[0], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[6]);
  FFV1P0_3(w[3], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[7]);
  // Amplitude(s) for diagram number 1
  VVVV1_0(w[6], w[7], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[0]); //printf(" &amp[0] %f %f\n", amp[0].real(), amp[0].imag());
  VVVV3_0(w[6], w[7], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[1]); //printf(" &amp[1] %f %f\n", amp[1].real(), amp[1].imag());
  VVVV4_0(w[6], w[7], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[2]); //printf(" &amp[2] %f %f\n", amp[2].real(), amp[2].imag());
  VVV1P0_1(w[6], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[8]);
  // Amplitude(s) for diagram number 2
  VVV1_0(w[7], w[5], w[8], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[3]); //printf(" &amp[3] %f %f\n", amp[3].real(), amp[3].imag());
  VVV1P0_1(w[6], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[9]);
  // Amplitude(s) for diagram number 3
  VVV1_0(w[7], w[4], w[9], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[4]); //printf(" &amp[4] %f %f\n", amp[4].real(), amp[4].imag());
  VVV1P0_1(w[4], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[10]);
  // Amplitude(s) for diagram number 4
  VVV1_0(w[6], w[7], w[10], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[5]); //printf(" &amp[5] %f %f\n", amp[5].real(), amp[5].imag());
  FFV1_1(w[2], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[11]);
  FFV1_2(w[3], w[6], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  // Amplitude(s) for diagram number 5
  FFV1_0(w[12], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[6]); //printf(" &amp[6] %f %f\n", amp[6].real(), amp[6].imag());
  // Amplitude(s) for diagram number 6
  FFV1_0(w[3], w[11], w[9], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[7]); //printf(" &amp[7] %f %f\n", amp[7].real(), amp[7].imag());
  FFV1_2(w[3], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[13]);
  // Amplitude(s) for diagram number 7
  FFV1_0(w[13], w[11], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[8]); //printf(" &amp[8] %f %f\n", amp[8].real(), amp[8].imag());
  FFV1_1(w[2], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[14]);
  // Amplitude(s) for diagram number 8
  FFV1_0(w[12], w[14], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[9]); //printf(" &amp[9] %f %f\n", amp[9].real(), amp[9].imag());
  // Amplitude(s) for diagram number 9
  FFV1_0(w[3], w[14], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[10]); //printf(" &amp[10] %f %f\n", amp[10].real(), amp[10].imag());
  FFV1_2(w[3], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[15]);
  // Amplitude(s) for diagram number 10
  FFV1_0(w[15], w[14], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[11]); //printf(" &amp[11] %f %f\n", amp[11].real(), amp[11].imag());
  FFV1_1(w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[16]);
  // Amplitude(s) for diagram number 11
  FFV1_0(w[15], w[16], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[12]); //printf(" &amp[12] %f %f\n", amp[12].real(), amp[12].imag());
  // Amplitude(s) for diagram number 12
  FFV1_0(w[15], w[2], w[9], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[13]); //printf(" &amp[13] %f %f\n", amp[13].real(), amp[13].imag());
  // Amplitude(s) for diagram number 13
  FFV1_0(w[13], w[16], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[14]); //printf(" &amp[14] %f %f\n", amp[14].real(), amp[14].imag());
  // Amplitude(s) for diagram number 14
  FFV1_0(w[13], w[2], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[15]); //printf(" &amp[15] %f %f\n", amp[15].real(), amp[15].imag());
  // Amplitude(s) for diagram number 15
  FFV1_0(w[3], w[16], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[16]); //printf(" &amp[16] %f %f\n", amp[16].real(), amp[16].imag());
  // Amplitude(s) for diagram number 16
  FFV1_0(w[12], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[17]); //printf(" &amp[17] %f %f\n", amp[17].real(), amp[17].imag());
  FFV1_1(w[2], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  FFV1_2(w[3], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[16]);
  FFV1_1(w[12], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[8]);
  // Amplitude(s) for diagram number 17
  FFV1_0(w[16], w[8], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[18]); //printf(" &amp[18] %f %f\n", amp[18].real(), amp[18].imag());
  FFV1_1(w[12], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[9]);
  // Amplitude(s) for diagram number 18
  FFV1_0(w[16], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[19]); //printf(" &amp[19] %f %f\n", amp[19].real(), amp[19].imag());
  // Amplitude(s) for diagram number 19
  FFV1_0(w[16], w[12], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[20]); //printf(" &amp[20] %f %f\n", amp[20].real(), amp[20].imag());
  VVV1P0_1(w[1], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[6]);
  FFV1P0_3(w[3], w[12], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[17]);
  // Amplitude(s) for diagram number 20
  VVV1_0(w[6], w[5], w[17], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[21]); //printf(" &amp[21] %f %f\n", amp[21].real(), amp[21].imag());
  // Amplitude(s) for diagram number 21
  FFV1_0(w[3], w[9], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[22]); //printf(" &amp[22] %f %f\n", amp[22].real(), amp[22].imag());
  // Amplitude(s) for diagram number 22
  FFV1_0(w[13], w[12], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[23]); //printf(" &amp[23] %f %f\n", amp[23].real(), amp[23].imag());
  VVV1P0_1(w[1], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[18]);
  // Amplitude(s) for diagram number 23
  VVV1_0(w[18], w[4], w[17], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[24]); //printf(" &amp[24] %f %f\n", amp[24].real(), amp[24].imag());
  // Amplitude(s) for diagram number 24
  FFV1_0(w[3], w[8], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[25]); //printf(" &amp[25] %f %f\n", amp[25].real(), amp[25].imag());
  // Amplitude(s) for diagram number 25
  FFV1_0(w[15], w[12], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[26]); //printf(" &amp[26] %f %f\n", amp[26].real(), amp[26].imag());
  FFV1_1(w[12], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[19]);
  // Amplitude(s) for diagram number 26
  FFV1_0(w[15], w[19], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[27]); //printf(" &amp[27] %f %f\n", amp[27].real(), amp[27].imag());
  // Amplitude(s) for diagram number 27
  FFV1_0(w[15], w[9], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[28]); //printf(" &amp[28] %f %f\n", amp[28].real(), amp[28].imag());
  // Amplitude(s) for diagram number 28
  FFV1_0(w[13], w[19], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[29]); //printf(" &amp[29] %f %f\n", amp[29].real(), amp[29].imag());
  // Amplitude(s) for diagram number 29
  FFV1_0(w[13], w[8], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[30]); //printf(" &amp[30] %f %f\n", amp[30].real(), amp[30].imag());
  // Amplitude(s) for diagram number 30
  FFV1_0(w[3], w[19], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[31]); //printf(" &amp[31] %f %f\n", amp[31].real(), amp[31].imag());
  // Amplitude(s) for diagram number 31
  VVV1_0(w[1], w[10], w[17], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[32]); //printf(" &amp[32] %f %f\n", amp[32].real(), amp[32].imag());
  VVVV1P0_1(w[1], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[17]);
  VVVV3P0_1(w[1], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[19]);
  VVVV4P0_1(w[1], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[8]);
  // Amplitude(s) for diagram number 32
  FFV1_0(w[3], w[12], w[17], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[33]); //printf(" &amp[33] %f %f\n", amp[33].real(), amp[33].imag());
  FFV1_0(w[3], w[12], w[19], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[34]); //printf(" &amp[34] %f %f\n", amp[34].real(), amp[34].imag());
  FFV1_0(w[3], w[12], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[35]); //printf(" &amp[35] %f %f\n", amp[35].real(), amp[35].imag());
  FFV1_2(w[3], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  FFV1_1(w[2], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[9]);
  FFV1_2(w[12], w[4], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[20]);
  // Amplitude(s) for diagram number 33
  FFV1_0(w[20], w[9], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[36]); //printf(" &amp[36] %f %f\n", amp[36].real(), amp[36].imag());
  FFV1_2(w[12], w[5], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[21]);
  // Amplitude(s) for diagram number 34
  FFV1_0(w[21], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[37]); //printf(" &amp[37] %f %f\n", amp[37].real(), amp[37].imag());
  // Amplitude(s) for diagram number 35
  FFV1_0(w[12], w[9], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[38]); //printf(" &amp[38] %f %f\n", amp[38].real(), amp[38].imag());
  FFV1P0_3(w[12], w[2], thrust::complex<double> (cIPC[2], cIPC[3]), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 36
  VVV1_0(w[6], w[5], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[39]); //printf(" &amp[39] %f %f\n", amp[39].real(), amp[39].imag());
  // Amplitude(s) for diagram number 37
  FFV1_0(w[21], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[40]); //printf(" &amp[40] %f %f\n", amp[40].real(), amp[40].imag());
  // Amplitude(s) for diagram number 38
  FFV1_0(w[12], w[14], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[41]); //printf(" &amp[41] %f %f\n", amp[41].real(), amp[41].imag());
  // Amplitude(s) for diagram number 39
  VVV1_0(w[18], w[4], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[42]); //printf(" &amp[42] %f %f\n", amp[42].real(), amp[42].imag());
  // Amplitude(s) for diagram number 40
  FFV1_0(w[20], w[2], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[43]); //printf(" &amp[43] %f %f\n", amp[43].real(), amp[43].imag());
  // Amplitude(s) for diagram number 41
  FFV1_0(w[12], w[11], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[44]); //printf(" &amp[44] %f %f\n", amp[44].real(), amp[44].imag());
  FFV1_2(w[12], w[1], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[23]);
  // Amplitude(s) for diagram number 42
  FFV1_0(w[23], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[45]); //printf(" &amp[45] %f %f\n", amp[45].real(), amp[45].imag());
  // Amplitude(s) for diagram number 43
  FFV1_0(w[21], w[11], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[46]); //printf(" &amp[46] %f %f\n", amp[46].real(), amp[46].imag());
  // Amplitude(s) for diagram number 44
  FFV1_0(w[23], w[14], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[47]); //printf(" &amp[47] %f %f\n", amp[47].real(), amp[47].imag());
  // Amplitude(s) for diagram number 45
  FFV1_0(w[20], w[14], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[48]); //printf(" &amp[48] %f %f\n", amp[48].real(), amp[48].imag());
  // Amplitude(s) for diagram number 46
  FFV1_0(w[23], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[49]); //printf(" &amp[49] %f %f\n", amp[49].real(), amp[49].imag());
  // Amplitude(s) for diagram number 47
  VVV1_0(w[1], w[10], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[50]); //printf(" &amp[50] %f %f\n", amp[50].real(), amp[50].imag());
  // Amplitude(s) for diagram number 48
  FFV1_0(w[12], w[2], w[17], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[51]); //printf(" &amp[51] %f %f\n", amp[51].real(), amp[51].imag());
  FFV1_0(w[12], w[2], w[19], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[52]); //printf(" &amp[52] %f %f\n", amp[52].real(), amp[52].imag());
  FFV1_0(w[12], w[2], w[8], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[53]); //printf(" &amp[53] %f %f\n", amp[53].real(), amp[53].imag());
  VVV1P0_1(w[0], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[12]);
  FFV1_2(w[3], w[12], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[22]);
  // Amplitude(s) for diagram number 49
  FFV1_0(w[22], w[9], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[54]); //printf(" &amp[54] %f %f\n", amp[54].real(), amp[54].imag());
  VVV1P0_1(w[12], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[23]);
  // Amplitude(s) for diagram number 50
  FFV1_0(w[3], w[9], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[55]); //printf(" &amp[55] %f %f\n", amp[55].real(), amp[55].imag());
  // Amplitude(s) for diagram number 51
  FFV1_0(w[13], w[9], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[56]); //printf(" &amp[56] %f %f\n", amp[56].real(), amp[56].imag());
  FFV1_1(w[2], w[12], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[20]);
  // Amplitude(s) for diagram number 52
  FFV1_0(w[16], w[20], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[57]); //printf(" &amp[57] %f %f\n", amp[57].real(), amp[57].imag());
  // Amplitude(s) for diagram number 53
  FFV1_0(w[16], w[2], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[58]); //printf(" &amp[58] %f %f\n", amp[58].real(), amp[58].imag());
  // Amplitude(s) for diagram number 54
  FFV1_0(w[16], w[14], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[59]); //printf(" &amp[59] %f %f\n", amp[59].real(), amp[59].imag());
  // Amplitude(s) for diagram number 55
  FFV1_0(w[3], w[20], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[60]); //printf(" &amp[60] %f %f\n", amp[60].real(), amp[60].imag());
  // Amplitude(s) for diagram number 56
  FFV1_0(w[22], w[2], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[61]); //printf(" &amp[61] %f %f\n", amp[61].real(), amp[61].imag());
  // Amplitude(s) for diagram number 57
  VVV1_0(w[12], w[18], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[62]); //printf(" &amp[62] %f %f\n", amp[62].real(), amp[62].imag());
  // Amplitude(s) for diagram number 58
  VVVV1_0(w[12], w[1], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[63]); //printf(" &amp[63] %f %f\n", amp[63].real(), amp[63].imag());
  VVVV3_0(w[12], w[1], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[64]); //printf(" &amp[64] %f %f\n", amp[64].real(), amp[64].imag());
  VVVV4_0(w[12], w[1], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[65]); //printf(" &amp[65] %f %f\n", amp[65].real(), amp[65].imag());
  VVV1P0_1(w[12], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[21]);
  // Amplitude(s) for diagram number 59
  VVV1_0(w[7], w[5], w[21], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[66]); //printf(" &amp[66] %f %f\n", amp[66].real(), amp[66].imag());
  // Amplitude(s) for diagram number 60
  VVV1_0(w[1], w[7], w[23], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[67]); //printf(" &amp[67] %f %f\n", amp[67].real(), amp[67].imag());
  // Amplitude(s) for diagram number 61
  FFV1_0(w[3], w[14], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[68]); //printf(" &amp[68] %f %f\n", amp[68].real(), amp[68].imag());
  // Amplitude(s) for diagram number 62
  FFV1_0(w[22], w[14], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[69]); //printf(" &amp[69] %f %f\n", amp[69].real(), amp[69].imag());
  // Amplitude(s) for diagram number 63
  FFV1_0(w[13], w[2], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[70]); //printf(" &amp[70] %f %f\n", amp[70].real(), amp[70].imag());
  // Amplitude(s) for diagram number 64
  FFV1_0(w[13], w[20], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[71]); //printf(" &amp[71] %f %f\n", amp[71].real(), amp[71].imag());
  VVV1P0_1(w[0], w[5], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[20]);
  FFV1_2(w[3], w[20], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[21]);
  // Amplitude(s) for diagram number 65
  FFV1_0(w[21], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[72]); //printf(" &amp[72] %f %f\n", amp[72].real(), amp[72].imag());
  VVV1P0_1(w[20], w[4], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 66
  FFV1_0(w[3], w[9], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[73]); //printf(" &amp[73] %f %f\n", amp[73].real(), amp[73].imag());
  // Amplitude(s) for diagram number 67
  FFV1_0(w[15], w[9], w[20], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[74]); //printf(" &amp[74] %f %f\n", amp[74].real(), amp[74].imag());
  FFV1_1(w[2], w[20], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[23]);
  // Amplitude(s) for diagram number 68
  FFV1_0(w[16], w[23], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[75]); //printf(" &amp[75] %f %f\n", amp[75].real(), amp[75].imag());
  // Amplitude(s) for diagram number 69
  FFV1_0(w[16], w[2], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[76]); //printf(" &amp[76] %f %f\n", amp[76].real(), amp[76].imag());
  // Amplitude(s) for diagram number 70
  FFV1_0(w[16], w[11], w[20], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[77]); //printf(" &amp[77] %f %f\n", amp[77].real(), amp[77].imag());
  // Amplitude(s) for diagram number 71
  FFV1_0(w[3], w[23], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[78]); //printf(" &amp[78] %f %f\n", amp[78].real(), amp[78].imag());
  // Amplitude(s) for diagram number 72
  FFV1_0(w[21], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[79]); //printf(" &amp[79] %f %f\n", amp[79].real(), amp[79].imag());
  // Amplitude(s) for diagram number 73
  VVV1_0(w[20], w[6], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[80]); //printf(" &amp[80] %f %f\n", amp[80].real(), amp[80].imag());
  // Amplitude(s) for diagram number 74
  VVVV1_0(w[20], w[1], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[81]); //printf(" &amp[81] %f %f\n", amp[81].real(), amp[81].imag());
  VVVV3_0(w[20], w[1], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[82]); //printf(" &amp[82] %f %f\n", amp[82].real(), amp[82].imag());
  VVVV4_0(w[20], w[1], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[83]); //printf(" &amp[83] %f %f\n", amp[83].real(), amp[83].imag());
  VVV1P0_1(w[20], w[1], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[12]);
  // Amplitude(s) for diagram number 75
  VVV1_0(w[7], w[4], w[12], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[84]); //printf(" &amp[84] %f %f\n", amp[84].real(), amp[84].imag());
  // Amplitude(s) for diagram number 76
  VVV1_0(w[1], w[7], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[85]); //printf(" &amp[85] %f %f\n", amp[85].real(), amp[85].imag());
  // Amplitude(s) for diagram number 77
  FFV1_0(w[3], w[11], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[86]); //printf(" &amp[86] %f %f\n", amp[86].real(), amp[86].imag());
  // Amplitude(s) for diagram number 78
  FFV1_0(w[21], w[11], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[87]); //printf(" &amp[87] %f %f\n", amp[87].real(), amp[87].imag());
  // Amplitude(s) for diagram number 79
  FFV1_0(w[15], w[2], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[88]); //printf(" &amp[88] %f %f\n", amp[88].real(), amp[88].imag());
  // Amplitude(s) for diagram number 80
  FFV1_0(w[15], w[23], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[89]); //printf(" &amp[89] %f %f\n", amp[89].real(), amp[89].imag());
  FFV1_1(w[9], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[23]);
  // Amplitude(s) for diagram number 81
  FFV1_0(w[15], w[23], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[90]); //printf(" &amp[90] %f %f\n", amp[90].real(), amp[90].imag());
  FFV1_2(w[15], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[12]);
  // Amplitude(s) for diagram number 82
  FFV1_0(w[12], w[9], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[91]); //printf(" &amp[91] %f %f\n", amp[91].real(), amp[91].imag());
  // Amplitude(s) for diagram number 83
  FFV1_0(w[13], w[23], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[92]); //printf(" &amp[92] %f %f\n", amp[92].real(), amp[92].imag());
  FFV1_2(w[13], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[21]);
  // Amplitude(s) for diagram number 84
  FFV1_0(w[21], w[9], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[93]); //printf(" &amp[93] %f %f\n", amp[93].real(), amp[93].imag());
  // Amplitude(s) for diagram number 85
  FFV1_0(w[3], w[23], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[94]); //printf(" &amp[94] %f %f\n", amp[94].real(), amp[94].imag());
  VVV1P0_1(w[0], w[10], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[23]);
  // Amplitude(s) for diagram number 86
  FFV1_0(w[3], w[9], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[95]); //printf(" &amp[95] %f %f\n", amp[95].real(), amp[95].imag());
  FFV1_2(w[16], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[22]);
  // Amplitude(s) for diagram number 87
  FFV1_0(w[22], w[11], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[96]); //printf(" &amp[96] %f %f\n", amp[96].real(), amp[96].imag());
  FFV1_1(w[11], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[20]);
  // Amplitude(s) for diagram number 88
  FFV1_0(w[16], w[20], w[5], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[97]); //printf(" &amp[97] %f %f\n", amp[97].real(), amp[97].imag());
  // Amplitude(s) for diagram number 89
  FFV1_0(w[22], w[14], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[98]); //printf(" &amp[98] %f %f\n", amp[98].real(), amp[98].imag());
  FFV1_1(w[14], w[0], thrust::complex<double> (cIPC[2], cIPC[3]), cIPD[0],
      cIPD[1], w[24]);
  // Amplitude(s) for diagram number 90
  FFV1_0(w[16], w[24], w[4], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[99]); //printf(" &amp[99] %f %f\n", amp[99].real(), amp[99].imag());
  // Amplitude(s) for diagram number 91
  FFV1_0(w[22], w[2], w[10], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[100]); //printf(" &amp[100] %f %f\n", amp[100].real(), amp[100].imag());
  // Amplitude(s) for diagram number 92
  FFV1_0(w[16], w[2], w[23], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[101]); //printf(" &amp[101] %f %f\n", amp[101].real(), amp[101].imag());
  // Amplitude(s) for diagram number 93
  VVVV1_0(w[0], w[6], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[102]); //printf(" &amp[102] %f %f\n", amp[102].real(), amp[102].imag());
  VVVV3_0(w[0], w[6], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[103]); //printf(" &amp[103] %f %f\n", amp[103].real(), amp[103].imag());
  VVVV4_0(w[0], w[6], w[7], w[5], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[104]); //printf(" &amp[104] %f %f\n", amp[104].real(), amp[104].imag());
  VVV1P0_1(w[0], w[6], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[22]);
  // Amplitude(s) for diagram number 94
  VVV1_0(w[7], w[5], w[22], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[105]); //printf(" &amp[105] %f %f\n", amp[105].real(), amp[105].imag());
  VVV1P0_1(w[0], w[7], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[25]);
  // Amplitude(s) for diagram number 95
  VVV1_0(w[6], w[5], w[25], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[106]); //printf(" &amp[106] %f %f\n", amp[106].real(), amp[106].imag());
  // Amplitude(s) for diagram number 96
  FFV1_0(w[3], w[14], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[107]); //printf(" &amp[107] %f %f\n", amp[107].real(), amp[107].imag());
  // Amplitude(s) for diagram number 97
  FFV1_0(w[3], w[24], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[108]); //printf(" &amp[108] %f %f\n", amp[108].real(), amp[108].imag());
  // Amplitude(s) for diagram number 98
  FFV1_0(w[13], w[2], w[22], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[109]); //printf(" &amp[109] %f %f\n", amp[109].real(), amp[109].imag());
  // Amplitude(s) for diagram number 99
  FFV1_0(w[21], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[110]); //printf(" &amp[110] %f %f\n", amp[110].real(), amp[110].imag());
  // Amplitude(s) for diagram number 100
  VVVV1_0(w[0], w[18], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[111]); //printf(" &amp[111] %f %f\n", amp[111].real(), amp[111].imag());
  VVVV3_0(w[0], w[18], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[112]); //printf(" &amp[112] %f %f\n", amp[112].real(), amp[112].imag());
  VVVV4_0(w[0], w[18], w[7], w[4], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[113]); //printf(" &amp[113] %f %f\n", amp[113].real(), amp[113].imag());
  VVV1P0_1(w[0], w[18], thrust::complex<double> (cIPC[0], cIPC[1]), 0., 0.,
      w[6]);
  // Amplitude(s) for diagram number 101
  VVV1_0(w[7], w[4], w[6], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[114]); //printf(" &amp[114] %f %f\n", amp[114].real(), amp[114].imag());
  // Amplitude(s) for diagram number 102
  VVV1_0(w[18], w[4], w[25], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[115]); //printf(" &amp[115] %f %f\n", amp[115].real(), amp[115].imag());
  // Amplitude(s) for diagram number 103
  FFV1_0(w[3], w[11], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[116]); //printf(" &amp[116] %f %f\n", amp[116].real(), amp[116].imag());
  // Amplitude(s) for diagram number 104
  FFV1_0(w[3], w[20], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[117]); //printf(" &amp[117] %f %f\n", amp[117].real(), amp[117].imag());
  // Amplitude(s) for diagram number 105
  FFV1_0(w[15], w[2], w[6], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[118]); //printf(" &amp[118] %f %f\n", amp[118].real(), amp[118].imag());
  // Amplitude(s) for diagram number 106
  FFV1_0(w[12], w[2], w[18], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[119]); //printf(" &amp[119] %f %f\n", amp[119].real(), amp[119].imag());
  // Amplitude(s) for diagram number 107
  VVVV1_0(w[0], w[1], w[7], w[10], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[120]); //printf(" &amp[120] %f %f\n", amp[120].real(), amp[120].imag());
  VVVV3_0(w[0], w[1], w[7], w[10], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[121]); //printf(" &amp[121] %f %f\n", amp[121].real(), amp[121].imag());
  VVVV4_0(w[0], w[1], w[7], w[10], thrust::complex<double> (cIPC[4], cIPC[5]),
      &amp[122]); //printf(" &amp[122] %f %f\n", amp[122].real(), amp[122].imag());
  // Amplitude(s) for diagram number 108
  VVV1_0(w[1], w[10], w[25], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[123]); //printf(" &amp[123] %f %f\n", amp[123].real(), amp[123].imag());
  // Amplitude(s) for diagram number 109
  VVV1_0(w[1], w[7], w[23], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[124]); //printf(" &amp[124] %f %f\n", amp[124].real(), amp[124].imag());
  // Amplitude(s) for diagram number 110
  FFV1_0(w[13], w[20], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[125]); //printf(" &amp[125] %f %f\n", amp[125].real(), amp[125].imag());
  // Amplitude(s) for diagram number 111
  FFV1_0(w[21], w[11], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[126]); //printf(" &amp[126] %f %f\n", amp[126].real(), amp[126].imag());
  // Amplitude(s) for diagram number 112
  FFV1_0(w[15], w[24], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[127]); //printf(" &amp[127] %f %f\n", amp[127].real(), amp[127].imag());
  // Amplitude(s) for diagram number 113
  FFV1_0(w[12], w[14], w[1], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[128]); //printf(" &amp[128] %f %f\n", amp[128].real(), amp[128].imag());
  VVVV1P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[12]);
  VVVV3P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[24]);
  VVVV4P0_1(w[0], w[1], w[4], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[21]);
  // Amplitude(s) for diagram number 114
  VVV1_0(w[12], w[7], w[5], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[129]); //printf(" &amp[129] %f %f\n", amp[129].real(), amp[129].imag());
  VVV1_0(w[24], w[7], w[5], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[130]); //printf(" &amp[130] %f %f\n", amp[130].real(), amp[130].imag());
  VVV1_0(w[21], w[7], w[5], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[131]); //printf(" &amp[131] %f %f\n", amp[131].real(), amp[131].imag());
  // Amplitude(s) for diagram number 115
  FFV1_0(w[3], w[14], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[132]); //printf(" &amp[132] %f %f\n", amp[132].real(), amp[132].imag());
  FFV1_0(w[3], w[14], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[133]); //printf(" &amp[133] %f %f\n", amp[133].real(), amp[133].imag());
  FFV1_0(w[3], w[14], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[134]); //printf(" &amp[134] %f %f\n", amp[134].real(), amp[134].imag());
  // Amplitude(s) for diagram number 116
  FFV1_0(w[13], w[2], w[12], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[135]); //printf(" &amp[135] %f %f\n", amp[135].real(), amp[135].imag());
  FFV1_0(w[13], w[2], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[136]); //printf(" &amp[136] %f %f\n", amp[136].real(), amp[136].imag());
  FFV1_0(w[13], w[2], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[137]); //printf(" &amp[137] %f %f\n", amp[137].real(), amp[137].imag());
  VVVV1P0_1(w[0], w[1], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[21]);
  VVVV3P0_1(w[0], w[1], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[13]);
  VVVV4P0_1(w[0], w[1], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[24]);
  // Amplitude(s) for diagram number 117
  VVV1_0(w[21], w[7], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[138]); //printf(" &amp[138] %f %f\n", amp[138].real(), amp[138].imag());
  VVV1_0(w[13], w[7], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[139]); //printf(" &amp[139] %f %f\n", amp[139].real(), amp[139].imag());
  VVV1_0(w[24], w[7], w[4], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[140]); //printf(" &amp[140] %f %f\n", amp[140].real(), amp[140].imag());
  // Amplitude(s) for diagram number 118
  FFV1_0(w[3], w[11], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[141]); //printf(" &amp[141] %f %f\n", amp[141].real(), amp[141].imag());
  FFV1_0(w[3], w[11], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[142]); //printf(" &amp[142] %f %f\n", amp[142].real(), amp[142].imag());
  FFV1_0(w[3], w[11], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[143]); //printf(" &amp[143] %f %f\n", amp[143].real(), amp[143].imag());
  // Amplitude(s) for diagram number 119
  FFV1_0(w[15], w[2], w[21], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[144]); //printf(" &amp[144] %f %f\n", amp[144].real(), amp[144].imag());
  FFV1_0(w[15], w[2], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[145]); //printf(" &amp[145] %f %f\n", amp[145].real(), amp[145].imag());
  FFV1_0(w[15], w[2], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[146]); //printf(" &amp[146] %f %f\n", amp[146].real(), amp[146].imag());
  VVVV1P0_1(w[0], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[24]);
  VVVV3P0_1(w[0], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[15]);
  VVVV4P0_1(w[0], w[4], w[5], thrust::complex<double> (cIPC[4], cIPC[5]), 0.,
      0., w[13]);
  // Amplitude(s) for diagram number 120
  FFV1_0(w[3], w[9], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[147]); //printf(" &amp[147] %f %f\n", amp[147].real(), amp[147].imag());
  FFV1_0(w[3], w[9], w[15], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[148]); //printf(" &amp[148] %f %f\n", amp[148].real(), amp[148].imag());
  FFV1_0(w[3], w[9], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[149]); //printf(" &amp[149] %f %f\n", amp[149].real(), amp[149].imag());
  // Amplitude(s) for diagram number 121
  FFV1_0(w[16], w[2], w[24], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[150]); //printf(" &amp[150] %f %f\n", amp[150].real(), amp[150].imag());
  FFV1_0(w[16], w[2], w[15], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[151]); //printf(" &amp[151] %f %f\n", amp[151].real(), amp[151].imag());
  FFV1_0(w[16], w[2], w[13], thrust::complex<double> (cIPC[2], cIPC[3]),
      &amp[152]); //printf(" &amp[152] %f %f\n", amp[152].real(), amp[152].imag());
  // Amplitude(s) for diagram number 122
  VVV1_0(w[24], w[1], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[153]); //printf(" &amp[153] %f %f\n", amp[153].real(), amp[153].imag());
  VVV1_0(w[15], w[1], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[154]); //printf(" &amp[154] %f %f\n", amp[154].real(), amp[154].imag());
  VVV1_0(w[13], w[1], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[155]); //printf(" &amp[155] %f %f\n", amp[155].real(), amp[155].imag());
  // Amplitude(s) for diagram number 123
  VVV1_0(w[0], w[17], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[156]); //printf(" &amp[156] %f %f\n", amp[156].real(), amp[156].imag());
  VVV1_0(w[0], w[19], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[157]); //printf(" &amp[157] %f %f\n", amp[157].real(), amp[157].imag());
  VVV1_0(w[0], w[8], w[7], thrust::complex<double> (cIPC[0], cIPC[1]),
      &amp[158]); //printf(" &amp[158] %f %f\n", amp[158].real(), amp[158].imag());
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
  size_t slicePitch = dpt * 6; 

  char * dps = devPtr + dim * slicePitch; 
  double * matrix_element = (double * )((char * )meDevPtr + dim * mePitch); 

  thrust::complex<double> amp[159]; 


  // Local variables and constants
  const int ncomb = 64; 
  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // std::complex<double> **wfs;
  double t[1]; 
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


