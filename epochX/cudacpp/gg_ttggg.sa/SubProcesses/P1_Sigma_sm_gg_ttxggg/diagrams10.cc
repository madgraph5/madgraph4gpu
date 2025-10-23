// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessChannelIds.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessWavefunctions.h"
#include "color_sum.h"
#include "diagrams.h"
#include "diagrams_header.h"

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
#include "MemoryAccessDenominators.h"
#include "MemoryAccessNumerators.h"
#endif

#ifdef MGONGPUCPP_GPUIMPL
namespace mg5amcGpu
#else
namespace mg5amcCpu
#endif
{
  //--------------------------------------------------------------------------

  __global__ void
  diagramgroup10( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                  fptype* jamps,                  // output jamps[ncolor*2*nevt]
                  const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                  cxtype_sv* jamp_sv,             // output jamps[ncolor*2*neppV]
                  const fptype** COUPs,           // input: dependent and independent COUPs[nxcoup] for this event page
#endif
                  const unsigned int* channelIds, // input: channelIds[nevt] for GPU or SCALAR channelId[0] for C++ (1 to #diagrams, 0 to disable SDE)
                  fptype* numerators,             // input/output: multichannel numerators[nevtORneppV], add helicity ihel
                  fptype* denominators,           // input/output: multichannel denominators[nevtORneppV], add helicity ihel
                  const fptype* cIPC,             // input: GPU __device__ or GPU host address of cIPC
                  const fptype* cIPD )            // input: GPU __device__ or GPU host address of cIPD
  {
    // A uniform interface for diagramgroupXXX including channelIDs, numerators and denominators is used also #ifndef MGONGPU_SUPPORTS_MULTICHANNEL
    // In that case, however, the boilerplate code asserts that all three pointers all nullptr as a sanity check
#include "diagrams_boilerplate.h"

#ifdef MGONGPUCPP_GPUIMPL
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) retrieveWf( wfs, w_cx, nevt, iwf );
    retrieveWf( wfs, w_cx, nevt, 0 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 23 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 25 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 33 );
    retrieveWf( wfs, w_cx, nevt, 35 );
    retrieveWf( wfs, w_cx, nevt, 36 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 41 );
    retrieveWf( wfs, w_cx, nevt, 43 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 49 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 117 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 119 );
    retrieveWf( wfs, w_cx, nevt, 120 );
#endif

    // *** DIAGRAM 901 OF 1240 ***
    // Wavefunction(s) for diagram number 901
    // (none)
    // Amplitude(s) for diagram number 901
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[69], w_fp[11], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 902 OF 1240 ***
    // Wavefunction(s) for diagram number 902
    // (none)
    // Amplitude(s) for diagram number 902
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[70], w_fp[13], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];

    // *** DIAGRAM 903 OF 1240 ***
    // Wavefunction(s) for diagram number 903
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[93] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[90] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 903
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[93], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[90], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 904 OF 1240 ***
    // Wavefunction(s) for diagram number 904
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[22] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[103] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[63] );
    // Amplitude(s) for diagram number 904
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[103], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[63], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 905 OF 1240 ***
    // Wavefunction(s) for diagram number 905
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[107] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[95] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[105] );
    // Amplitude(s) for diagram number 905
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[6], w_fp[107], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[6], w_fp[95], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[6], w_fp[105], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 906 OF 1240 ***
    // Wavefunction(s) for diagram number 906
    // (none)
    // Amplitude(s) for diagram number 906
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[118], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[119], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[120], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 907 OF 1240 ***
    // Wavefunction(s) for diagram number 907
    // (none)
    // Amplitude(s) for diagram number 907
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[8], w_fp[27], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[8], w_fp[27], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], w_fp[8], w_fp[27], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 908 OF 1240 ***
    // Wavefunction(s) for diagram number 908
    // (none)
    // Amplitude(s) for diagram number 908
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[27], w_fp[65], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 909 OF 1240 ***
    // Wavefunction(s) for diagram number 909
    // (none)
    // Amplitude(s) for diagram number 909
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[27], w_fp[56], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 910 OF 1240 ***
    // Wavefunction(s) for diagram number 910
    // (none)
    // Amplitude(s) for diagram number 910
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[8], w_fp[101], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 911 OF 1240 ***
    // Wavefunction(s) for diagram number 911
    // (none)
    // Amplitude(s) for diagram number 911
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[37], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 912 OF 1240 ***
    // Wavefunction(s) for diagram number 912
    // (none)
    // Amplitude(s) for diagram number 912
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[36], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];

    // *** DIAGRAM 913 OF 1240 ***
    // Wavefunction(s) for diagram number 913
    // (none)
    // Amplitude(s) for diagram number 913
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[114], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 914 OF 1240 ***
    // Wavefunction(s) for diagram number 914
    // (none)
    // Amplitude(s) for diagram number 914
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[114], w_fp[70], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];

    // *** DIAGRAM 915 OF 1240 ***
    // Wavefunction(s) for diagram number 915
    // (none)
    // Amplitude(s) for diagram number 915
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[36], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 916 OF 1240 ***
    // Wavefunction(s) for diagram number 916
    // (none)
    // Amplitude(s) for diagram number 916
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[70], w_fp[37], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 917 OF 1240 ***
    // Wavefunction(s) for diagram number 917
    // (none)
    // Amplitude(s) for diagram number 917
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 918 OF 1240 ***
    // Wavefunction(s) for diagram number 918
    // (none)
    // Amplitude(s) for diagram number 918
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[33], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];

    // *** DIAGRAM 919 OF 1240 ***
    // Wavefunction(s) for diagram number 919
    // (none)
    // Amplitude(s) for diagram number 919
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[114], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 920 OF 1240 ***
    // Wavefunction(s) for diagram number 920
    // (none)
    // Amplitude(s) for diagram number 920
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[33], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 921 OF 1240 ***
    // Wavefunction(s) for diagram number 921
    // (none)
    // Amplitude(s) for diagram number 921
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[51], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 922 OF 1240 ***
    // Wavefunction(s) for diagram number 922
    // (none)
    // Amplitude(s) for diagram number 922
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[49], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 923 OF 1240 ***
    // Wavefunction(s) for diagram number 923
    // (none)
    // Amplitude(s) for diagram number 923
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[113], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 924 OF 1240 ***
    // Wavefunction(s) for diagram number 924
    // (none)
    // Amplitude(s) for diagram number 924
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[113], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];

    // *** DIAGRAM 925 OF 1240 ***
    // Wavefunction(s) for diagram number 925
    // (none)
    // Amplitude(s) for diagram number 925
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[49], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 926 OF 1240 ***
    // Wavefunction(s) for diagram number 926
    // (none)
    // Amplitude(s) for diagram number 926
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[69], w_fp[51], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 927 OF 1240 ***
    // Wavefunction(s) for diagram number 927
    // (none)
    // Amplitude(s) for diagram number 927
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[93], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 928 OF 1240 ***
    // Wavefunction(s) for diagram number 928
    // (none)
    // Amplitude(s) for diagram number 928
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[47], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 929 OF 1240 ***
    // Wavefunction(s) for diagram number 929
    // (none)
    // Amplitude(s) for diagram number 929
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[113], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 930 OF 1240 ***
    // Wavefunction(s) for diagram number 930
    // (none)
    // Amplitude(s) for diagram number 930
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[47], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 931 OF 1240 ***
    // Wavefunction(s) for diagram number 931
    // (none)
    // Amplitude(s) for diagram number 931
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[54], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 932 OF 1240 ***
    // Wavefunction(s) for diagram number 932
    // (none)
    // Amplitude(s) for diagram number 932
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[2], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];

    // *** DIAGRAM 933 OF 1240 ***
    // Wavefunction(s) for diagram number 933
    // (none)
    // Amplitude(s) for diagram number 933
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[94], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 934 OF 1240 ***
    // Wavefunction(s) for diagram number 934
    // (none)
    // Amplitude(s) for diagram number 934
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[2], w_fp[70], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 935 OF 1240 ***
    // Wavefunction(s) for diagram number 935
    // (none)
    // Amplitude(s) for diagram number 935
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[94], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 936 OF 1240 ***
    // Wavefunction(s) for diagram number 936
    // (none)
    // Amplitude(s) for diagram number 936
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[70], w_fp[54], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 937 OF 1240 ***
    // Wavefunction(s) for diagram number 937
    // (none)
    // Amplitude(s) for diagram number 937
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 938 OF 1240 ***
    // Wavefunction(s) for diagram number 938
    // (none)
    // Amplitude(s) for diagram number 938
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[20], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 939 OF 1240 ***
    // Wavefunction(s) for diagram number 939
    // (none)
    // Amplitude(s) for diagram number 939
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[28], w_fp[2], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];

    // *** DIAGRAM 940 OF 1240 ***
    // Wavefunction(s) for diagram number 940
    // (none)
    // Amplitude(s) for diagram number 940
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[94], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 941 OF 1240 ***
    // Wavefunction(s) for diagram number 941
    // (none)
    // Amplitude(s) for diagram number 941
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[2], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];

    // *** DIAGRAM 942 OF 1240 ***
    // Wavefunction(s) for diagram number 942
    // (none)
    // Amplitude(s) for diagram number 942
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[28], w_fp[94], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 943 OF 1240 ***
    // Wavefunction(s) for diagram number 943
    // (none)
    // Amplitude(s) for diagram number 943
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[69], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 944 OF 1240 ***
    // Wavefunction(s) for diagram number 944
    // (none)
    // Amplitude(s) for diagram number 944
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[93], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 945 OF 1240 ***
    // Wavefunction(s) for diagram number 945
    // (none)
    // Amplitude(s) for diagram number 945
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[15], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 946 OF 1240 ***
    // Wavefunction(s) for diagram number 946
    // (none)
    // Amplitude(s) for diagram number 946
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[2], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 947 OF 1240 ***
    // Wavefunction(s) for diagram number 947
    // (none)
    // Amplitude(s) for diagram number 947
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[94], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 948 OF 1240 ***
    // Wavefunction(s) for diagram number 948
    // (none)
    // Amplitude(s) for diagram number 948
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[2], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 949 OF 1240 ***
    // Wavefunction(s) for diagram number 949
    // (none)
    // Amplitude(s) for diagram number 949
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[94], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];

    // *** DIAGRAM 950 OF 1240 ***
    // Wavefunction(s) for diagram number 950
    // (none)
    // Amplitude(s) for diagram number 950
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[15], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 951 OF 1240 ***
    // Wavefunction(s) for diagram number 951
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], COUPs[0], 1.0, 0., 0., w_fp[71] );
    // Amplitude(s) for diagram number 951
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[13], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 952 OF 1240 ***
    // Wavefunction(s) for diagram number 952
    // (none)
    // Amplitude(s) for diagram number 952
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[10], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];

    // *** DIAGRAM 953 OF 1240 ***
    // Wavefunction(s) for diagram number 953
    // (none)
    // Amplitude(s) for diagram number 953
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[5], w_fp[71], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[5], w_fp[71], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[5], w_fp[71], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 954 OF 1240 ***
    // Wavefunction(s) for diagram number 954
    // (none)
    // Amplitude(s) for diagram number 954
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[74], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];

    // *** DIAGRAM 955 OF 1240 ***
    // Wavefunction(s) for diagram number 955
    // (none)
    // Amplitude(s) for diagram number 955
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[75], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 956 OF 1240 ***
    // Wavefunction(s) for diagram number 956
    // (none)
    // Amplitude(s) for diagram number 956
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 957 OF 1240 ***
    // Wavefunction(s) for diagram number 957
    // (none)
    // Amplitude(s) for diagram number 957
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[74], w_fp[10], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];

    // *** DIAGRAM 958 OF 1240 ***
    // Wavefunction(s) for diagram number 958
    // (none)
    // Amplitude(s) for diagram number 958
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[75], w_fp[13], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];

    // *** DIAGRAM 959 OF 1240 ***
    // Wavefunction(s) for diagram number 959
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[94] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[65] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 959
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[94], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[65], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];

    // *** DIAGRAM 960 OF 1240 ***
    // Wavefunction(s) for diagram number 960
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[90] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[93] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[69] );
    // Amplitude(s) for diagram number 960
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[90], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[93], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[69], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 961 OF 1240 ***
    // Wavefunction(s) for diagram number 961
    // (none)
    // Amplitude(s) for diagram number 961
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[5], w_fp[107], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[5], w_fp[95], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[5], w_fp[105], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 962 OF 1240 ***
    // Wavefunction(s) for diagram number 962
    // (none)
    // Amplitude(s) for diagram number 962
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[115], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[117], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];

    // *** DIAGRAM 963 OF 1240 ***
    // Wavefunction(s) for diagram number 963
    // (none)
    // Amplitude(s) for diagram number 963
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[8], w_fp[24], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[8], w_fp[24], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[72], w_fp[8], w_fp[24], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 964 OF 1240 ***
    // Wavefunction(s) for diagram number 964
    // (none)
    // Amplitude(s) for diagram number 964
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[24], w_fp[71], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 965 OF 1240 ***
    // Wavefunction(s) for diagram number 965
    // (none)
    // Amplitude(s) for diagram number 965
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[24], w_fp[56], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 966 OF 1240 ***
    // Wavefunction(s) for diagram number 966
    // (none)
    // Amplitude(s) for diagram number 966
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[8], w_fp[98], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 967 OF 1240 ***
    // Wavefunction(s) for diagram number 967
    // (none)
    // Amplitude(s) for diagram number 967
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[37], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 968 OF 1240 ***
    // Wavefunction(s) for diagram number 968
    // (none)
    // Amplitude(s) for diagram number 968
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[35], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];

    // *** DIAGRAM 969 OF 1240 ***
    // Wavefunction(s) for diagram number 969
    // (none)
    // Amplitude(s) for diagram number 969
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[76], w_fp[114], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 970 OF 1240 ***
    // Wavefunction(s) for diagram number 970
    // (none)
    // Amplitude(s) for diagram number 970
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[114], w_fp[75], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];

    // *** DIAGRAM 971 OF 1240 ***
    // Wavefunction(s) for diagram number 971
    // (none)
    // Amplitude(s) for diagram number 971
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[76], w_fp[35], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 972 OF 1240 ***
    // Wavefunction(s) for diagram number 972
    // (none)
    // Amplitude(s) for diagram number 972
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[75], w_fp[37], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 973 OF 1240 ***
    // Wavefunction(s) for diagram number 973
    // (none)
    // Amplitude(s) for diagram number 973
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[93], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 974 OF 1240 ***
    // Wavefunction(s) for diagram number 974
    // (none)
    // Amplitude(s) for diagram number 974
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[33], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];

    // *** DIAGRAM 975 OF 1240 ***
    // Wavefunction(s) for diagram number 975
    // (none)
    // Amplitude(s) for diagram number 975
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[114], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 976 OF 1240 ***
    // Wavefunction(s) for diagram number 976
    // (none)
    // Amplitude(s) for diagram number 976
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[33], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 977 OF 1240 ***
    // Wavefunction(s) for diagram number 977
    // (none)
    // Amplitude(s) for diagram number 977
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[45], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 978 OF 1240 ***
    // Wavefunction(s) for diagram number 978
    // (none)
    // Amplitude(s) for diagram number 978
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[43], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];

    // *** DIAGRAM 979 OF 1240 ***
    // Wavefunction(s) for diagram number 979
    // (none)
    // Amplitude(s) for diagram number 979
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[76], w_fp[102], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 980 OF 1240 ***
    // Wavefunction(s) for diagram number 980
    // (none)
    // Amplitude(s) for diagram number 980
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[102], w_fp[74], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];

    // *** DIAGRAM 981 OF 1240 ***
    // Wavefunction(s) for diagram number 981
    // (none)
    // Amplitude(s) for diagram number 981
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[76], w_fp[43], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 982 OF 1240 ***
    // Wavefunction(s) for diagram number 982
    // (none)
    // Amplitude(s) for diagram number 982
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[74], w_fp[45], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 983 OF 1240 ***
    // Wavefunction(s) for diagram number 983
    // (none)
    // Amplitude(s) for diagram number 983
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[94], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 984 OF 1240 ***
    // Wavefunction(s) for diagram number 984
    // (none)
    // Amplitude(s) for diagram number 984
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[39], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];

    // *** DIAGRAM 985 OF 1240 ***
    // Wavefunction(s) for diagram number 985
    // (none)
    // Amplitude(s) for diagram number 985
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[102], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 986 OF 1240 ***
    // Wavefunction(s) for diagram number 986
    // (none)
    // Amplitude(s) for diagram number 986
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[39], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 987 OF 1240 ***
    // Wavefunction(s) for diagram number 987
    // (none)
    // Amplitude(s) for diagram number 987
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[54], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 988 OF 1240 ***
    // Wavefunction(s) for diagram number 988
    // (none)
    // Amplitude(s) for diagram number 988
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[2], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 989 OF 1240 ***
    // Wavefunction(s) for diagram number 989
    // (none)
    // Amplitude(s) for diagram number 989
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[97], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 990 OF 1240 ***
    // Wavefunction(s) for diagram number 990
    // (none)
    // Amplitude(s) for diagram number 990
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[2], w_fp[75], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];

    // *** DIAGRAM 991 OF 1240 ***
    // Wavefunction(s) for diagram number 991
    // (none)
    // Amplitude(s) for diagram number 991
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[97], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 992 OF 1240 ***
    // Wavefunction(s) for diagram number 992
    // (none)
    // Amplitude(s) for diagram number 992
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[75], w_fp[54], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 993 OF 1240 ***
    // Wavefunction(s) for diagram number 993
    // (none)
    // Amplitude(s) for diagram number 993
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[93], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 994 OF 1240 ***
    // Wavefunction(s) for diagram number 994
    // (none)
    // Amplitude(s) for diagram number 994
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[23], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 995 OF 1240 ***
    // Wavefunction(s) for diagram number 995
    // (none)
    // Amplitude(s) for diagram number 995
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[25], w_fp[2], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];

    // *** DIAGRAM 996 OF 1240 ***
    // Wavefunction(s) for diagram number 996
    // (none)
    // Amplitude(s) for diagram number 996
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[97], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 997 OF 1240 ***
    // Wavefunction(s) for diagram number 997
    // (none)
    // Amplitude(s) for diagram number 997
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[2], w_fp[74], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];

    // *** DIAGRAM 998 OF 1240 ***
    // Wavefunction(s) for diagram number 998
    // (none)
    // Amplitude(s) for diagram number 998
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[25], w_fp[97], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 999 OF 1240 ***
    // Wavefunction(s) for diagram number 999
    // (none)
    // Amplitude(s) for diagram number 999
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[74], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1000 OF 1240 ***
    // Wavefunction(s) for diagram number 1000
    // (none)
    // Amplitude(s) for diagram number 1000
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[94], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[65], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 21 );
    storeWf( wfs, w_cx, nevt, 22 );
    storeWf( wfs, w_cx, nevt, 63 );
    storeWf( wfs, w_cx, nevt, 65 );
    storeWf( wfs, w_cx, nevt, 69 );
    storeWf( wfs, w_cx, nevt, 71 );
    storeWf( wfs, w_cx, nevt, 90 );
    storeWf( wfs, w_cx, nevt, 93 );
    storeWf( wfs, w_cx, nevt, 94 );
    storeWf( wfs, w_cx, nevt, 95 );
    storeWf( wfs, w_cx, nevt, 103 );
    storeWf( wfs, w_cx, nevt, 105 );
    storeWf( wfs, w_cx, nevt, 107 );
#endif
  }

  //--------------------------------------------------------------------------
}
