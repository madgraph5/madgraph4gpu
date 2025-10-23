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
  diagramgroup12( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 23 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 33 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 41 );
    retrieveWf( wfs, w_cx, nevt, 42 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 78 );
    retrieveWf( wfs, w_cx, nevt, 96 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 111 );
#endif

    // *** DIAGRAM 1101 OF 1240 ***
    // Wavefunction(s) for diagram number 1101
    // (none)
    // Amplitude(s) for diagram number 1101
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[59], w_fp[51], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1101 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1102 OF 1240 ***
    // Wavefunction(s) for diagram number 1102
    // (none)
    // Amplitude(s) for diagram number 1102
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[97], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1103 OF 1240 ***
    // Wavefunction(s) for diagram number 1103
    // (none)
    // Amplitude(s) for diagram number 1103
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[2], w_fp[67], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1103 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1104 OF 1240 ***
    // Wavefunction(s) for diagram number 1104
    // (none)
    // Amplitude(s) for diagram number 1104
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[18], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1104 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1105 OF 1240 ***
    // Wavefunction(s) for diagram number 1105
    // (none)
    // Amplitude(s) for diagram number 1105
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[78], w_fp[2], w_fp[96], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1105 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];

    // *** DIAGRAM 1106 OF 1240 ***
    // Wavefunction(s) for diagram number 1106
    // (none)
    // Amplitude(s) for diagram number 1106
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[96], w_fp[1], w_fp[54], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1106 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1107 OF 1240 ***
    // Wavefunction(s) for diagram number 1107
    // (none)
    // Amplitude(s) for diagram number 1107
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[78], w_fp[18], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1107 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1108 OF 1240 ***
    // Wavefunction(s) for diagram number 1108
    // (none)
    // Amplitude(s) for diagram number 1108
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[67], w_fp[54], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1108 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1109 OF 1240 ***
    // Wavefunction(s) for diagram number 1109
    // (none)
    // Amplitude(s) for diagram number 1109
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[76], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[42], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[11], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1110 OF 1240 ***
    // Wavefunction(s) for diagram number 1110
    // (none)
    // Amplitude(s) for diagram number 1110
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1110 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 1111 OF 1240 ***
    // Wavefunction(s) for diagram number 1111
    // (none)
    // Amplitude(s) for diagram number 1111
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[15], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1111 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1112 OF 1240 ***
    // Wavefunction(s) for diagram number 1112
    // (none)
    // Amplitude(s) for diagram number 1112
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[58], w_fp[2], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1112 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];

    // *** DIAGRAM 1113 OF 1240 ***
    // Wavefunction(s) for diagram number 1113
    // (none)
    // Amplitude(s) for diagram number 1113
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[101], w_fp[1], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1113 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1114 OF 1240 ***
    // Wavefunction(s) for diagram number 1114
    // (none)
    // Amplitude(s) for diagram number 1114
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[58], w_fp[15], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1114 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1115 OF 1240 ***
    // Wavefunction(s) for diagram number 1115
    // (none)
    // Amplitude(s) for diagram number 1115
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[68], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1115 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1116 OF 1240 ***
    // Wavefunction(s) for diagram number 1116
    // (none)
    // Amplitude(s) for diagram number 1116
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[111], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1117 OF 1240 ***
    // Wavefunction(s) for diagram number 1117
    // (none)
    // Amplitude(s) for diagram number 1117
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[2], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1117 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 1118 OF 1240 ***
    // Wavefunction(s) for diagram number 1118
    // (none)
    // Amplitude(s) for diagram number 1118
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[17], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1118 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1119 OF 1240 ***
    // Wavefunction(s) for diagram number 1119
    // (none)
    // Amplitude(s) for diagram number 1119
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[2], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1119 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];

    // *** DIAGRAM 1120 OF 1240 ***
    // Wavefunction(s) for diagram number 1120
    // (none)
    // Amplitude(s) for diagram number 1120
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[1], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1120 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1121 OF 1240 ***
    // Wavefunction(s) for diagram number 1121
    // (none)
    // Amplitude(s) for diagram number 1121
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[17], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1121 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1122 OF 1240 ***
    // Wavefunction(s) for diagram number 1122
    // (none)
    // Amplitude(s) for diagram number 1122
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[59], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 1122 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1123 OF 1240 ***
    // Wavefunction(s) for diagram number 1123
    // (none)
    // Amplitude(s) for diagram number 1123
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[97], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1124 OF 1240 ***
    // Wavefunction(s) for diagram number 1124
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[21] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[71] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[4], COUPs[2], 1.0, 0., 0., w_fp[97] );
    // Amplitude(s) for diagram number 1124
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[8], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 1125 OF 1240 ***
    // Wavefunction(s) for diagram number 1125
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[21], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[59] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[71], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[20] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[97], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[60] );
    // Amplitude(s) for diagram number 1125
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[59], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[60], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1126 OF 1240 ***
    // Wavefunction(s) for diagram number 1126
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[21], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[17] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[71], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[98] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[97], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[111] );
    // Amplitude(s) for diagram number 1126
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[98], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[111], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 1127 OF 1240 ***
    // Wavefunction(s) for diagram number 1127
    // (none)
    // Amplitude(s) for diagram number 1127
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[8], w_fp[29], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[8], w_fp[29], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[8], w_fp[29], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1128 OF 1240 ***
    // Wavefunction(s) for diagram number 1128
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[21], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[71], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[10] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[97], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[68] );
    // Amplitude(s) for diagram number 1128
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[39], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[39], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[39], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];

    // *** DIAGRAM 1129 OF 1240 ***
    // Wavefunction(s) for diagram number 1129
    // (none)
    // Amplitude(s) for diagram number 1129
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[111], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1130 OF 1240 ***
    // Wavefunction(s) for diagram number 1130
    // (none)
    // Amplitude(s) for diagram number 1130
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[39], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[39], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[39], w_fp[97], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];

    // *** DIAGRAM 1131 OF 1240 ***
    // Wavefunction(s) for diagram number 1131
    // (none)
    // Amplitude(s) for diagram number 1131
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[47], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[47], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[47], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1132 OF 1240 ***
    // Wavefunction(s) for diagram number 1132
    // (none)
    // Amplitude(s) for diagram number 1132
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1133 OF 1240 ***
    // Wavefunction(s) for diagram number 1133
    // (none)
    // Amplitude(s) for diagram number 1133
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[47], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[47], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[47], w_fp[97], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];

    // *** DIAGRAM 1134 OF 1240 ***
    // Wavefunction(s) for diagram number 1134
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[21], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[71], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[97], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[71] );
    // Amplitude(s) for diagram number 1134
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[23], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[21], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[71], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];

    // *** DIAGRAM 1135 OF 1240 ***
    // Wavefunction(s) for diagram number 1135
    // (none)
    // Amplitude(s) for diagram number 1135
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[111], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1136 OF 1240 ***
    // Wavefunction(s) for diagram number 1136
    // (none)
    // Amplitude(s) for diagram number 1136
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[23], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[21], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[71], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];

    // *** DIAGRAM 1137 OF 1240 ***
    // Wavefunction(s) for diagram number 1137
    // (none)
    // Amplitude(s) for diagram number 1137
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1138 OF 1240 ***
    // Wavefunction(s) for diagram number 1138
    // (none)
    // Amplitude(s) for diagram number 1138
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[21], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[71], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1139 OF 1240 ***
    // Wavefunction(s) for diagram number 1139
    // (none)
    // Amplitude(s) for diagram number 1139
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1140 OF 1240 ***
    // Wavefunction(s) for diagram number 1140
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[68] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[29] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 1140
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[8], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1141 OF 1240 ***
    // Wavefunction(s) for diagram number 1141
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[68], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[16] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[29], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[71] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[10], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 1141
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[71], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1142 OF 1240 ***
    // Wavefunction(s) for diagram number 1142
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[68], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[23] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[29], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[60] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[10], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[20] );
    // Amplitude(s) for diagram number 1142
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[60], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1143 OF 1240 ***
    // Wavefunction(s) for diagram number 1143
    // (none)
    // Amplitude(s) for diagram number 1143
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[68], w_fp[8], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
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
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[29], w_fp[8], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[10], w_fp[8], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 1144 OF 1240 ***
    // Wavefunction(s) for diagram number 1144
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[68], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[59] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[29], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[111] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[10], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[98] );
    // Amplitude(s) for diagram number 1144
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[33], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[33], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[33], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];

    // *** DIAGRAM 1145 OF 1240 ***
    // Wavefunction(s) for diagram number 1145
    // (none)
    // Amplitude(s) for diagram number 1145
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1146 OF 1240 ***
    // Wavefunction(s) for diagram number 1146
    // (none)
    // Amplitude(s) for diagram number 1146
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[33], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[33], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[33], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];

    // *** DIAGRAM 1147 OF 1240 ***
    // Wavefunction(s) for diagram number 1147
    // (none)
    // Amplitude(s) for diagram number 1147
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[47], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[47], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[47], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 1148 OF 1240 ***
    // Wavefunction(s) for diagram number 1148
    // (none)
    // Amplitude(s) for diagram number 1148
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1149 OF 1240 ***
    // Wavefunction(s) for diagram number 1149
    // (none)
    // Amplitude(s) for diagram number 1149
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[47], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[47], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[47], w_fp[10], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];

    // *** DIAGRAM 1150 OF 1240 ***
    // Wavefunction(s) for diagram number 1150
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[68], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[17] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[29], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[68] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[10], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[29] );
    // Amplitude(s) for diagram number 1150
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[17], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[68], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[29], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];

    // *** DIAGRAM 1151 OF 1240 ***
    // Wavefunction(s) for diagram number 1151
    // (none)
    // Amplitude(s) for diagram number 1151
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1152 OF 1240 ***
    // Wavefunction(s) for diagram number 1152
    // (none)
    // Amplitude(s) for diagram number 1152
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[17], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[68], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[29], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];

    // *** DIAGRAM 1153 OF 1240 ***
    // Wavefunction(s) for diagram number 1153
    // (none)
    // Amplitude(s) for diagram number 1153
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1154 OF 1240 ***
    // Wavefunction(s) for diagram number 1154
    // (none)
    // Amplitude(s) for diagram number 1154
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[17], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[68], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[29], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1155 OF 1240 ***
    // Wavefunction(s) for diagram number 1155
    // (none)
    // Amplitude(s) for diagram number 1155
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[2], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[2], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[2], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1156 OF 1240 ***
    // Wavefunction(s) for diagram number 1156
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[98] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[27] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[111] );
    // Amplitude(s) for diagram number 1156
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[8], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 1157 OF 1240 ***
    // Wavefunction(s) for diagram number 1157
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[98], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[59] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[27], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[29] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[111], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[68] );
    // Amplitude(s) for diagram number 1157
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[59], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[29], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[68], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];

    // *** DIAGRAM 1158 OF 1240 ***
    // Wavefunction(s) for diagram number 1158
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[98], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[17] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[27], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[21] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[111], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[71] );
    // Amplitude(s) for diagram number 1158
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[17], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[71], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 1159 OF 1240 ***
    // Wavefunction(s) for diagram number 1159
    // (none)
    // Amplitude(s) for diagram number 1159
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[8], w_fp[24], COUPs[0], 1.0, &amp_fp[0] );
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
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[8], w_fp[24], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[111], w_fp[8], w_fp[24], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 1160 OF 1240 ***
    // Wavefunction(s) for diagram number 1160
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[27], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[20] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[111], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[60] );
    // Amplitude(s) for diagram number 1160
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[33], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[33], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[33], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];

    // *** DIAGRAM 1161 OF 1240 ***
    // Wavefunction(s) for diagram number 1161
    // (none)
    // Amplitude(s) for diagram number 1161
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[33], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1162 OF 1240 ***
    // Wavefunction(s) for diagram number 1162
    // (none)
    // Amplitude(s) for diagram number 1162
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[33], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[33], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[33], w_fp[111], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];

    // *** DIAGRAM 1163 OF 1240 ***
    // Wavefunction(s) for diagram number 1163
    // (none)
    // Amplitude(s) for diagram number 1163
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[39], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[39], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[39], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];

    // *** DIAGRAM 1164 OF 1240 ***
    // Wavefunction(s) for diagram number 1164
    // (none)
    // Amplitude(s) for diagram number 1164
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1165 OF 1240 ***
    // Wavefunction(s) for diagram number 1165
    // (none)
    // Amplitude(s) for diagram number 1165
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[39], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[39], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[39], w_fp[111], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];

    // *** DIAGRAM 1166 OF 1240 ***
    // Wavefunction(s) for diagram number 1166
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[98], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[27], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[98] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[111], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[27] );
    // Amplitude(s) for diagram number 1166
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[23], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[98], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[27], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 1167 OF 1240 ***
    // Wavefunction(s) for diagram number 1167
    // (none)
    // Amplitude(s) for diagram number 1167
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[2], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1168 OF 1240 ***
    // Wavefunction(s) for diagram number 1168
    // (none)
    // Amplitude(s) for diagram number 1168
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[23], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[98], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[27], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];

    // *** DIAGRAM 1169 OF 1240 ***
    // Wavefunction(s) for diagram number 1169
    // (none)
    // Amplitude(s) for diagram number 1169
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1170 OF 1240 ***
    // Wavefunction(s) for diagram number 1170
    // (none)
    // Amplitude(s) for diagram number 1170
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[27], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1171 OF 1240 ***
    // Wavefunction(s) for diagram number 1171
    // (none)
    // Amplitude(s) for diagram number 1171
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1172 OF 1240 ***
    // Wavefunction(s) for diagram number 1172
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[60] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[24] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[20] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[60], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[16] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[24], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[27] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[20], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[98] );
    // Amplitude(s) for diagram number 1172
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];

    // *** DIAGRAM 1173 OF 1240 ***
    // Wavefunction(s) for diagram number 1173
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[60], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[23] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[24], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[68] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[29] );
    // Amplitude(s) for diagram number 1173
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1174 OF 1240 ***
    // Wavefunction(s) for diagram number 1174
    // (none)
    // Amplitude(s) for diagram number 1174
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[77], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[77], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[77], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];

    // *** DIAGRAM 1175 OF 1240 ***
    // Wavefunction(s) for diagram number 1175
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[60], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[59] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[24], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[71] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[20], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 1175
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[59], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[71], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[21], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];

    // *** DIAGRAM 1176 OF 1240 ***
    // Wavefunction(s) for diagram number 1176
    // (none)
    // Amplitude(s) for diagram number 1176
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[23], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1177 OF 1240 ***
    // Wavefunction(s) for diagram number 1177
    // (none)
    // Amplitude(s) for diagram number 1177
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[47], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[47], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[47], w_fp[20], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1178 OF 1240 ***
    // Wavefunction(s) for diagram number 1178
    // (none)
    // Amplitude(s) for diagram number 1178
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[59], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[71], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[21], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1179 OF 1240 ***
    // Wavefunction(s) for diagram number 1179
    // (none)
    // Amplitude(s) for diagram number 1179
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[2], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[2], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[2], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1180 OF 1240 ***
    // Wavefunction(s) for diagram number 1180
    // (none)
    // Amplitude(s) for diagram number 1180
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[72], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
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
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[72], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[72], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 1181 OF 1240 ***
    // Wavefunction(s) for diagram number 1181
    // (none)
    // Amplitude(s) for diagram number 1181
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[15] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], w_fp[8], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[9] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1182 OF 1240 ***
    // Wavefunction(s) for diagram number 1182
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[60], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[72] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[24], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[60] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[20], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[24] );
    // Amplitude(s) for diagram number 1182
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[72], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[60], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[24], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1183 OF 1240 ***
    // Wavefunction(s) for diagram number 1183
    // (none)
    // Amplitude(s) for diagram number 1183
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[68], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[15] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[29], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[9] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1184 OF 1240 ***
    // Wavefunction(s) for diagram number 1184
    // (none)
    // Amplitude(s) for diagram number 1184
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1185 OF 1240 ***
    // Wavefunction(s) for diagram number 1185
    // (none)
    // Amplitude(s) for diagram number 1185
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[47], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[47], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[47], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 1186 OF 1240 ***
    // Wavefunction(s) for diagram number 1186
    // (none)
    // Amplitude(s) for diagram number 1186
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[60], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1187 OF 1240 ***
    // Wavefunction(s) for diagram number 1187
    // (none)
    // Amplitude(s) for diagram number 1187
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[59], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[71], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[21], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];

    // *** DIAGRAM 1188 OF 1240 ***
    // Wavefunction(s) for diagram number 1188
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[21] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[71] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[59] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[21], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[24] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[71], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[60] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[59], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[72] );
    // Amplitude(s) for diagram number 1188
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];

    // *** DIAGRAM 1189 OF 1240 ***
    // Wavefunction(s) for diagram number 1189
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[21], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[98] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[71], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[27] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[59], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[16] );
    // Amplitude(s) for diagram number 1189
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1190 OF 1240 ***
    // Wavefunction(s) for diagram number 1190
    // (none)
    // Amplitude(s) for diagram number 1190
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[77], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[77], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[77], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];

    // *** DIAGRAM 1191 OF 1240 ***
    // Wavefunction(s) for diagram number 1191
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[21], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[29] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[71], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[68] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[59], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[23] );
    // Amplitude(s) for diagram number 1191
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[29], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[68], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[23], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];

    // *** DIAGRAM 1192 OF 1240 ***
    // Wavefunction(s) for diagram number 1192
    // (none)
    // Amplitude(s) for diagram number 1192
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1193 OF 1240 ***
    // Wavefunction(s) for diagram number 1193
    // (none)
    // Amplitude(s) for diagram number 1193
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[39], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[39], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[39], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];

    // *** DIAGRAM 1194 OF 1240 ***
    // Wavefunction(s) for diagram number 1194
    // (none)
    // Amplitude(s) for diagram number 1194
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[29], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[68], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[23], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1195 OF 1240 ***
    // Wavefunction(s) for diagram number 1195
    // (none)
    // Amplitude(s) for diagram number 1195
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[2], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[2], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[2], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1196 OF 1240 ***
    // Wavefunction(s) for diagram number 1196
    // (none)
    // Amplitude(s) for diagram number 1196
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[66], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
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
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[66], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[66], w_fp[8], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];

    // *** DIAGRAM 1197 OF 1240 ***
    // Wavefunction(s) for diagram number 1197
    // (none)
    // Amplitude(s) for diagram number 1197
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[21] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[59], w_fp[1], w_fp[8], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[11] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 1198 OF 1240 ***
    // Wavefunction(s) for diagram number 1198
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[21], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[66] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[71], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[21] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[59], w_fp[1], COUPs[0], 1.0, 0., 0., w_fp[71] );
    // Amplitude(s) for diagram number 1198
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[66], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[71], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];

    // *** DIAGRAM 1199 OF 1240 ***
    // Wavefunction(s) for diagram number 1199
    // (none)
    // Amplitude(s) for diagram number 1199
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[98], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[27], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[21] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[8], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[11] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 1200 OF 1240 ***
    // Wavefunction(s) for diagram number 1200
    // (none)
    // Amplitude(s) for diagram number 1200
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[71], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 10 );
    storeWf( wfs, w_cx, nevt, 16 );
    storeWf( wfs, w_cx, nevt, 17 );
    storeWf( wfs, w_cx, nevt, 20 );
    storeWf( wfs, w_cx, nevt, 21 );
    storeWf( wfs, w_cx, nevt, 23 );
    storeWf( wfs, w_cx, nevt, 24 );
    storeWf( wfs, w_cx, nevt, 27 );
    storeWf( wfs, w_cx, nevt, 29 );
    storeWf( wfs, w_cx, nevt, 59 );
    storeWf( wfs, w_cx, nevt, 60 );
    storeWf( wfs, w_cx, nevt, 66 );
    storeWf( wfs, w_cx, nevt, 68 );
    storeWf( wfs, w_cx, nevt, 71 );
    storeWf( wfs, w_cx, nevt, 72 );
    storeWf( wfs, w_cx, nevt, 97 );
    storeWf( wfs, w_cx, nevt, 98 );
    storeWf( wfs, w_cx, nevt, 111 );
#endif
  }

  //--------------------------------------------------------------------------
}
