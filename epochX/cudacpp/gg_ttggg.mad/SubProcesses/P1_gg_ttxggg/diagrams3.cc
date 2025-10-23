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
  diagramgroup3( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 25 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 31 );
    retrieveWf( wfs, w_cx, nevt, 32 );
    retrieveWf( wfs, w_cx, nevt, 33 );
    retrieveWf( wfs, w_cx, nevt, 34 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 41 );
    retrieveWf( wfs, w_cx, nevt, 42 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 77 );
#endif

    // *** DIAGRAM 201 OF 1240 ***
    // Wavefunction(s) for diagram number 201
    // (none)
    // Amplitude(s) for diagram number 201
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[58], w_fp[57], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 201 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[20] -= amp_sv[0];

    // *** DIAGRAM 202 OF 1240 ***
    // Wavefunction(s) for diagram number 202
    // (none)
    // Amplitude(s) for diagram number 202
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[25], w_fp[57], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 202 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] -= amp_sv[0];

    // *** DIAGRAM 203 OF 1240 ***
    // Wavefunction(s) for diagram number 203
    // (none)
    // Amplitude(s) for diagram number 203
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[77], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 203 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 204 OF 1240 ***
    // Wavefunction(s) for diagram number 204
    // (none)
    // Amplitude(s) for diagram number 204
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], w_fp[67], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 204 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];

    // *** DIAGRAM 205 OF 1240 ***
    // Wavefunction(s) for diagram number 205
    // (none)
    // Amplitude(s) for diagram number 205
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[60], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 205 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 206 OF 1240 ***
    // Wavefunction(s) for diagram number 206
    // (none)
    // Amplitude(s) for diagram number 206
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[28], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 206 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= amp_sv[0];

    // *** DIAGRAM 207 OF 1240 ***
    // Wavefunction(s) for diagram number 207
    // (none)
    // Amplitude(s) for diagram number 207
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[77], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 207 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= amp_sv[0];

    // *** DIAGRAM 208 OF 1240 ***
    // Wavefunction(s) for diagram number 208
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[41], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[60] );
    // Amplitude(s) for diagram number 208
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 208 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] -= amp_sv[0];

    // *** DIAGRAM 209 OF 1240 ***
    // Wavefunction(s) for diagram number 209
    // (none)
    // Amplitude(s) for diagram number 209
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[9], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 209 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] -= amp_sv[0];

    // *** DIAGRAM 210 OF 1240 ***
    // Wavefunction(s) for diagram number 210
    // (none)
    // Amplitude(s) for diagram number 210
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[60], w_fp[55], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 210 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[14] -= amp_sv[0];

    // *** DIAGRAM 211 OF 1240 ***
    // Wavefunction(s) for diagram number 211
    // (none)
    // Amplitude(s) for diagram number 211
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[28], w_fp[55], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 211 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] -= amp_sv[0];

    // *** DIAGRAM 212 OF 1240 ***
    // Wavefunction(s) for diagram number 212
    // (none)
    // Amplitude(s) for diagram number 212
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[77], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 212 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 213 OF 1240 ***
    // Wavefunction(s) for diagram number 213
    // (none)
    // Amplitude(s) for diagram number 213
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[24], w_fp[68], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 213 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];

    // *** DIAGRAM 214 OF 1240 ***
    // Wavefunction(s) for diagram number 214
    // (none)
    // Amplitude(s) for diagram number 214
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[59], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 214 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 215 OF 1240 ***
    // Wavefunction(s) for diagram number 215
    // (none)
    // Amplitude(s) for diagram number 215
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[26], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 215 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 216 OF 1240 ***
    // Wavefunction(s) for diagram number 216
    // (none)
    // Amplitude(s) for diagram number 216
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[42], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 216 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];

    // *** DIAGRAM 217 OF 1240 ***
    // Wavefunction(s) for diagram number 217
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[24], COUPs[0], 1.0, 0., 0., w_fp[59] );
    // Amplitude(s) for diagram number 217
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[59], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 217 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 218 OF 1240 ***
    // Wavefunction(s) for diagram number 218
    // (none)
    // Amplitude(s) for diagram number 218
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[1], w_fp[42], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 218 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 219 OF 1240 ***
    // Wavefunction(s) for diagram number 219
    // (none)
    // Amplitude(s) for diagram number 219
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[24], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[24], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[24], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 220 OF 1240 ***
    // Wavefunction(s) for diagram number 220
    // (none)
    // Amplitude(s) for diagram number 220
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[57], w_fp[59], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 220 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 221 OF 1240 ***
    // Wavefunction(s) for diagram number 221
    // (none)
    // Amplitude(s) for diagram number 221
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[26], w_fp[57], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 221 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 222 OF 1240 ***
    // Wavefunction(s) for diagram number 222
    // (none)
    // Amplitude(s) for diagram number 222
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 222 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 223 OF 1240 ***
    // Wavefunction(s) for diagram number 223
    // (none)
    // Amplitude(s) for diagram number 223
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 223 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];

    // *** DIAGRAM 224 OF 1240 ***
    // Wavefunction(s) for diagram number 224
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], COUPs[0], 1.0, 0., 0., w_fp[68] );
    // Amplitude(s) for diagram number 224
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[68], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 224 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 225 OF 1240 ***
    // Wavefunction(s) for diagram number 225
    // (none)
    // Amplitude(s) for diagram number 225
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[1], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 225 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 226 OF 1240 ***
    // Wavefunction(s) for diagram number 226
    // (none)
    // Amplitude(s) for diagram number 226
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], w_fp[5], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], w_fp[5], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[27], w_fp[5], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 227 OF 1240 ***
    // Wavefunction(s) for diagram number 227
    // (none)
    // Amplitude(s) for diagram number 227
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[55], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 227 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];

    // *** DIAGRAM 228 OF 1240 ***
    // Wavefunction(s) for diagram number 228
    // (none)
    // Amplitude(s) for diagram number 228
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[55], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 228 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 229 OF 1240 ***
    // Wavefunction(s) for diagram number 229
    // (none)
    // Amplitude(s) for diagram number 229
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[77], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 229 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 230 OF 1240 ***
    // Wavefunction(s) for diagram number 230
    // (none)
    // Amplitude(s) for diagram number 230
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 230 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];

    // *** DIAGRAM 231 OF 1240 ***
    // Wavefunction(s) for diagram number 231
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[29], COUPs[0], 1.0, 0., 0., w_fp[67] );
    // Amplitude(s) for diagram number 231
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[67], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 231 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 232 OF 1240 ***
    // Wavefunction(s) for diagram number 232
    // (none)
    // Amplitude(s) for diagram number 232
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[1], w_fp[19], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 232 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 233 OF 1240 ***
    // Wavefunction(s) for diagram number 233
    // (none)
    // Amplitude(s) for diagram number 233
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[29], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[29], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[29], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 234 OF 1240 ***
    // Wavefunction(s) for diagram number 234
    // (none)
    // Amplitude(s) for diagram number 234
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[67], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 234 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];

    // *** DIAGRAM 235 OF 1240 ***
    // Wavefunction(s) for diagram number 235
    // (none)
    // Amplitude(s) for diagram number 235
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[9], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 235 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 236 OF 1240 ***
    // Wavefunction(s) for diagram number 236
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[73] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[79] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[80] );
    // Amplitude(s) for diagram number 236
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[73], w_fp[6], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[79], w_fp[6], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[80], w_fp[6], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 237 OF 1240 ***
    // Wavefunction(s) for diagram number 237
    // (none)
    // Amplitude(s) for diagram number 237
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[57], w_fp[73], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[57], w_fp[79], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[57], w_fp[80], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 238 OF 1240 ***
    // Wavefunction(s) for diagram number 238
    // (none)
    // Amplitude(s) for diagram number 238
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[34], w_fp[73], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[34], w_fp[79], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[34], w_fp[80], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];

    // *** DIAGRAM 239 OF 1240 ***
    // Wavefunction(s) for diagram number 239
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[57] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[81] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[82] );
    // Amplitude(s) for diagram number 239
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[57], w_fp[5], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[81], w_fp[5], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[82], w_fp[5], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 240 OF 1240 ***
    // Wavefunction(s) for diagram number 240
    // (none)
    // Amplitude(s) for diagram number 240
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[55], w_fp[57], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[55], w_fp[81], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[55], w_fp[82], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];

    // *** DIAGRAM 241 OF 1240 ***
    // Wavefunction(s) for diagram number 241
    // (none)
    // Amplitude(s) for diagram number 241
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[34], w_fp[57], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[34], w_fp[81], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[34], w_fp[82], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];

    // *** DIAGRAM 242 OF 1240 ***
    // Wavefunction(s) for diagram number 242
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[55] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[83] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[84] );
    // Amplitude(s) for diagram number 242
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[55], w_fp[4], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[83], w_fp[4], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[84], w_fp[4], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 243 OF 1240 ***
    // Wavefunction(s) for diagram number 243
    // (none)
    // Amplitude(s) for diagram number 243
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[55], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[83], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];

    // *** DIAGRAM 244 OF 1240 ***
    // Wavefunction(s) for diagram number 244
    // (none)
    // Amplitude(s) for diagram number 244
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[34], w_fp[55], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[34], w_fp[83], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[34], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];

    // *** DIAGRAM 245 OF 1240 ***
    // Wavefunction(s) for diagram number 245
    // (none)
    // Amplitude(s) for diagram number 245
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[30], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[31], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[32], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];

    // *** DIAGRAM 246 OF 1240 ***
    // Wavefunction(s) for diagram number 246
    // (none)
    // Amplitude(s) for diagram number 246
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[30], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[31], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[32], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 247 OF 1240 ***
    // Wavefunction(s) for diagram number 247
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[62] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[77] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[34] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[77], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    // Amplitude(s) for diagram number 247
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[9], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 247 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[41] -= amp_sv[0];

    // *** DIAGRAM 248 OF 1240 ***
    // Wavefunction(s) for diagram number 248
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[77], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[85] );
    // Amplitude(s) for diagram number 248
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[85], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 248 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[47] -= amp_sv[0];

    // *** DIAGRAM 249 OF 1240 ***
    // Wavefunction(s) for diagram number 249
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[86] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[77], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[87] );
    // Amplitude(s) for diagram number 249
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[87], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 249 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[35] -= amp_sv[0];

    // *** DIAGRAM 250 OF 1240 ***
    // Wavefunction(s) for diagram number 250
    // (none)
    // Amplitude(s) for diagram number 250
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[85], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 250 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[45] -= amp_sv[0];

    // *** DIAGRAM 251 OF 1240 ***
    // Wavefunction(s) for diagram number 251
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[88] );
    // Amplitude(s) for diagram number 251
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[87], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 251 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] -= amp_sv[0];

    // *** DIAGRAM 252 OF 1240 ***
    // Wavefunction(s) for diagram number 252
    // (none)
    // Amplitude(s) for diagram number 252
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[9], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 252 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] -= amp_sv[0];

    // *** DIAGRAM 253 OF 1240 ***
    // Wavefunction(s) for diagram number 253
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[77], COUPs[1], 1.0, 0., 0., w_fp[89] );
    // Amplitude(s) for diagram number 253
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[6], w_fp[89], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 253 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];

    // *** DIAGRAM 254 OF 1240 ***
    // Wavefunction(s) for diagram number 254
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[24], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[90] );
    // Amplitude(s) for diagram number 254
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 254 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 255 OF 1240 ***
    // Wavefunction(s) for diagram number 255
    // (none)
    // Amplitude(s) for diagram number 255
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[77], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 255 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 256 OF 1240 ***
    // Wavefunction(s) for diagram number 256
    // (none)
    // Amplitude(s) for diagram number 256
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[5], w_fp[89], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 256 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];

    // *** DIAGRAM 257 OF 1240 ***
    // Wavefunction(s) for diagram number 257
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[27], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[91] );
    // Amplitude(s) for diagram number 257
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[91], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 257 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 258 OF 1240 ***
    // Wavefunction(s) for diagram number 258
    // (none)
    // Amplitude(s) for diagram number 258
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[77], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 258 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 259 OF 1240 ***
    // Wavefunction(s) for diagram number 259
    // (none)
    // Amplitude(s) for diagram number 259
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[29], w_fp[89], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 259 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];

    // *** DIAGRAM 260 OF 1240 ***
    // Wavefunction(s) for diagram number 260
    // (none)
    // Amplitude(s) for diagram number 260
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[77], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 260 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 261 OF 1240 ***
    // Wavefunction(s) for diagram number 261
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[29], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[89] );
    // Amplitude(s) for diagram number 261
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[89], w_fp[77], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 261 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 262 OF 1240 ***
    // Wavefunction(s) for diagram number 262
    // (none)
    // Amplitude(s) for diagram number 262
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[77], w_fp[30], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[77], w_fp[31], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[77], w_fp[32], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];

    // *** DIAGRAM 263 OF 1240 ***
    // Wavefunction(s) for diagram number 263
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[2], COUPs[1], 1.0, 0., 0., w_fp[92] );
    // Amplitude(s) for diagram number 263
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[63], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 263 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 264 OF 1240 ***
    // Wavefunction(s) for diagram number 264
    // (none)
    // Amplitude(s) for diagram number 264
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[64], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 264 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 265 OF 1240 ***
    // Wavefunction(s) for diagram number 265
    // (none)
    // Amplitude(s) for diagram number 265
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 266 OF 1240 ***
    // Wavefunction(s) for diagram number 266
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[61], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[93] );
    // Amplitude(s) for diagram number 266
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[93], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 266 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 267 OF 1240 ***
    // Wavefunction(s) for diagram number 267
    // (none)
    // Amplitude(s) for diagram number 267
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[2], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 267 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];

    // *** DIAGRAM 268 OF 1240 ***
    // Wavefunction(s) for diagram number 268
    // (none)
    // Amplitude(s) for diagram number 268
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[93], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 268 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 269 OF 1240 ***
    // Wavefunction(s) for diagram number 269
    // (none)
    // Amplitude(s) for diagram number 269
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[2], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 269 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];

    // *** DIAGRAM 270 OF 1240 ***
    // Wavefunction(s) for diagram number 270
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[61], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[94] );
    // Amplitude(s) for diagram number 270
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[39], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 270 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 271 OF 1240 ***
    // Wavefunction(s) for diagram number 271
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[39], COUPs[1], 1.0, 0., 0., w_fp[95] );
    // Amplitude(s) for diagram number 271
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[6], w_fp[95], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 271 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];

    // *** DIAGRAM 272 OF 1240 ***
    // Wavefunction(s) for diagram number 272
    // (none)
    // Amplitude(s) for diagram number 272
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[39], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 272 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 273 OF 1240 ***
    // Wavefunction(s) for diagram number 273
    // (none)
    // Amplitude(s) for diagram number 273
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[47], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 273 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 274 OF 1240 ***
    // Wavefunction(s) for diagram number 274
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[47], COUPs[1], 1.0, 0., 0., w_fp[96] );
    // Amplitude(s) for diagram number 274
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[96], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 274 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 275 OF 1240 ***
    // Wavefunction(s) for diagram number 275
    // (none)
    // Amplitude(s) for diagram number 275
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[47], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 275 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 276 OF 1240 ***
    // Wavefunction(s) for diagram number 276
    // (none)
    // Amplitude(s) for diagram number 276
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[94], w_fp[2], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 276 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 277 OF 1240 ***
    // Wavefunction(s) for diagram number 277
    // (none)
    // Amplitude(s) for diagram number 277
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[29], w_fp[92], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 277 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 278 OF 1240 ***
    // Wavefunction(s) for diagram number 278
    // (none)
    // Amplitude(s) for diagram number 278
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[89], w_fp[2], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 278 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];

    // *** DIAGRAM 279 OF 1240 ***
    // Wavefunction(s) for diagram number 279
    // (none)
    // Amplitude(s) for diagram number 279
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[69], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 279 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 280 OF 1240 ***
    // Wavefunction(s) for diagram number 280
    // (none)
    // Amplitude(s) for diagram number 280
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[70], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 280 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 281 OF 1240 ***
    // Wavefunction(s) for diagram number 281
    // (none)
    // Amplitude(s) for diagram number 281
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 282 OF 1240 ***
    // Wavefunction(s) for diagram number 282
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[66], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[94] );
    // Amplitude(s) for diagram number 282
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[94], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 282 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 283 OF 1240 ***
    // Wavefunction(s) for diagram number 283
    // (none)
    // Amplitude(s) for diagram number 283
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[70], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 283 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];

    // *** DIAGRAM 284 OF 1240 ***
    // Wavefunction(s) for diagram number 284
    // (none)
    // Amplitude(s) for diagram number 284
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[94], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 284 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 285 OF 1240 ***
    // Wavefunction(s) for diagram number 285
    // (none)
    // Amplitude(s) for diagram number 285
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[2], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 285 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];

    // *** DIAGRAM 286 OF 1240 ***
    // Wavefunction(s) for diagram number 286
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[66], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[97] );
    // Amplitude(s) for diagram number 286
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[33], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 286 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 287 OF 1240 ***
    // Wavefunction(s) for diagram number 287
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[62], w_fp[33], COUPs[1], 1.0, 0., 0., w_fp[98] );
    // Amplitude(s) for diagram number 287
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[6], w_fp[98], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 287 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];

    // *** DIAGRAM 288 OF 1240 ***
    // Wavefunction(s) for diagram number 288
    // (none)
    // Amplitude(s) for diagram number 288
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[88], w_fp[33], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 288 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 289 OF 1240 ***
    // Wavefunction(s) for diagram number 289
    // (none)
    // Amplitude(s) for diagram number 289
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[47], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 289 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 290 OF 1240 ***
    // Wavefunction(s) for diagram number 290
    // (none)
    // Amplitude(s) for diagram number 290
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[96], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 290 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 291 OF 1240 ***
    // Wavefunction(s) for diagram number 291
    // (none)
    // Amplitude(s) for diagram number 291
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[47], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 291 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 292 OF 1240 ***
    // Wavefunction(s) for diagram number 292
    // (none)
    // Amplitude(s) for diagram number 292
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[97], w_fp[2], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 292 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 293 OF 1240 ***
    // Wavefunction(s) for diagram number 293
    // (none)
    // Amplitude(s) for diagram number 293
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[27], w_fp[92], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 293 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 294 OF 1240 ***
    // Wavefunction(s) for diagram number 294
    // (none)
    // Amplitude(s) for diagram number 294
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[91], w_fp[2], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 294 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];

    // *** DIAGRAM 295 OF 1240 ***
    // Wavefunction(s) for diagram number 295
    // (none)
    // Amplitude(s) for diagram number 295
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[74], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 295 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 296 OF 1240 ***
    // Wavefunction(s) for diagram number 296
    // (none)
    // Amplitude(s) for diagram number 296
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[92], w_fp[75], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 296 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 297 OF 1240 ***
    // Wavefunction(s) for diagram number 297
    // (none)
    // Amplitude(s) for diagram number 297
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[92], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 298 OF 1240 ***
    // Wavefunction(s) for diagram number 298
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[72], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[97] );
    // Amplitude(s) for diagram number 298
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[97], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 298 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 299 OF 1240 ***
    // Wavefunction(s) for diagram number 299
    // (none)
    // Amplitude(s) for diagram number 299
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[75], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 299 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[47] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 300 OF 1240 ***
    // Wavefunction(s) for diagram number 300
    // (none)
    // Amplitude(s) for diagram number 300
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[86], w_fp[97], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 300 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 34 );
    storeWf( wfs, w_cx, nevt, 55 );
    storeWf( wfs, w_cx, nevt, 57 );
    storeWf( wfs, w_cx, nevt, 59 );
    storeWf( wfs, w_cx, nevt, 60 );
    storeWf( wfs, w_cx, nevt, 62 );
    storeWf( wfs, w_cx, nevt, 67 );
    storeWf( wfs, w_cx, nevt, 68 );
    storeWf( wfs, w_cx, nevt, 73 );
    storeWf( wfs, w_cx, nevt, 77 );
    storeWf( wfs, w_cx, nevt, 79 );
    storeWf( wfs, w_cx, nevt, 80 );
    storeWf( wfs, w_cx, nevt, 81 );
    storeWf( wfs, w_cx, nevt, 82 );
    storeWf( wfs, w_cx, nevt, 83 );
    storeWf( wfs, w_cx, nevt, 84 );
    storeWf( wfs, w_cx, nevt, 85 );
    storeWf( wfs, w_cx, nevt, 86 );
    storeWf( wfs, w_cx, nevt, 87 );
    storeWf( wfs, w_cx, nevt, 88 );
    storeWf( wfs, w_cx, nevt, 89 );
    storeWf( wfs, w_cx, nevt, 90 );
    storeWf( wfs, w_cx, nevt, 91 );
    storeWf( wfs, w_cx, nevt, 92 );
    storeWf( wfs, w_cx, nevt, 93 );
    storeWf( wfs, w_cx, nevt, 94 );
    storeWf( wfs, w_cx, nevt, 95 );
    storeWf( wfs, w_cx, nevt, 96 );
    storeWf( wfs, w_cx, nevt, 97 );
    storeWf( wfs, w_cx, nevt, 98 );
#endif
  }

  //--------------------------------------------------------------------------
}
