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
  diagramgroup2( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 25 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 31 );
    retrieveWf( wfs, w_cx, nevt, 32 );
    retrieveWf( wfs, w_cx, nevt, 34 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 41 );
    retrieveWf( wfs, w_cx, nevt, 42 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 48 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
#endif

    // *** DIAGRAM 101 OF 1240 ***
    // Wavefunction(s) for diagram number 101
    // (none)
    // Amplitude(s) for diagram number 101
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[52], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 101 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];

    // *** DIAGRAM 102 OF 1240 ***
    // Wavefunction(s) for diagram number 102
    // (none)
    // Amplitude(s) for diagram number 102
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[42], w_fp[2], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 102 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 103 OF 1240 ***
    // Wavefunction(s) for diagram number 103
    // (none)
    // Amplitude(s) for diagram number 103
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[26], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 103 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 104 OF 1240 ***
    // Wavefunction(s) for diagram number 104
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[24], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[26] );
    // Amplitude(s) for diagram number 104
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[26], w_fp[52], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 104 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];

    // *** DIAGRAM 105 OF 1240 ***
    // Wavefunction(s) for diagram number 105
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[24], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[42] );
    // Amplitude(s) for diagram number 105
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[42], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 105 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 106 OF 1240 ***
    // Wavefunction(s) for diagram number 106
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[24], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[17] );
    // Amplitude(s) for diagram number 106
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[17], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 106 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];

    // *** DIAGRAM 107 OF 1240 ***
    // Wavefunction(s) for diagram number 107
    // (none)
    // Amplitude(s) for diagram number 107
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[42], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 107 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 108 OF 1240 ***
    // Wavefunction(s) for diagram number 108
    // (none)
    // Amplitude(s) for diagram number 108
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[17], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 108 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 109 OF 1240 ***
    // Wavefunction(s) for diagram number 109
    // (none)
    // Amplitude(s) for diagram number 109
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[26], w_fp[2], w_fp[14], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 109 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 110 OF 1240 ***
    // Wavefunction(s) for diagram number 110
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[27], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[14] );
    // Amplitude(s) for diagram number 110
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[52], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 110 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];

    // *** DIAGRAM 111 OF 1240 ***
    // Wavefunction(s) for diagram number 111
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[27], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[16] );
    // Amplitude(s) for diagram number 111
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 111 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 112 OF 1240 ***
    // Wavefunction(s) for diagram number 112
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[27], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[15] );
    // Amplitude(s) for diagram number 112
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[15], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 112 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 113 OF 1240 ***
    // Wavefunction(s) for diagram number 113
    // (none)
    // Amplitude(s) for diagram number 113
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 113 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 114 OF 1240 ***
    // Wavefunction(s) for diagram number 114
    // (none)
    // Amplitude(s) for diagram number 114
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[15], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 114 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 115 OF 1240 ***
    // Wavefunction(s) for diagram number 115
    // (none)
    // Amplitude(s) for diagram number 115
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[14], w_fp[2], w_fp[12], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 115 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 116 OF 1240 ***
    // Wavefunction(s) for diagram number 116
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[29], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 116
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[52], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 116 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];

    // *** DIAGRAM 117 OF 1240 ***
    // Wavefunction(s) for diagram number 117
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[4], w_fp[29], COUPs[0], 1.0, 0., 0., w_fp[19] );
    // Amplitude(s) for diagram number 117
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 117 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 118 OF 1240 ***
    // Wavefunction(s) for diagram number 118
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[29], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[18] );
    // Amplitude(s) for diagram number 118
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[18], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 118 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 119 OF 1240 ***
    // Wavefunction(s) for diagram number 119
    // (none)
    // Amplitude(s) for diagram number 119
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 119 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 120 OF 1240 ***
    // Wavefunction(s) for diagram number 120
    // (none)
    // Amplitude(s) for diagram number 120
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[18], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 120 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 121 OF 1240 ***
    // Wavefunction(s) for diagram number 121
    // (none)
    // Amplitude(s) for diagram number 121
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 121 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 122 OF 1240 ***
    // Wavefunction(s) for diagram number 122
    // (none)
    // Amplitude(s) for diagram number 122
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[30], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[31], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[52], w_fp[32], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 123 OF 1240 ***
    // Wavefunction(s) for diagram number 123
    // (none)
    // Amplitude(s) for diagram number 123
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[30], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[31], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[34], w_fp[2], w_fp[32], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 124 OF 1240 ***
    // Wavefunction(s) for diagram number 124
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[2], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[34] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[52] );
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[9] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[52], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[22] );
    // Amplitude(s) for diagram number 124
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[9], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 124 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] -= amp_sv[0];

    // *** DIAGRAM 125 OF 1240 ***
    // Wavefunction(s) for diagram number 125
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[52], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[21] );
    // Amplitude(s) for diagram number 125
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 125 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] -= amp_sv[0];

    // *** DIAGRAM 126 OF 1240 ***
    // Wavefunction(s) for diagram number 126
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[5], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[55] );
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[52], w_fp[4], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[56] );
    // Amplitude(s) for diagram number 126
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[55], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 126 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[17] -= amp_sv[0];

    // *** DIAGRAM 127 OF 1240 ***
    // Wavefunction(s) for diagram number 127
    // (none)
    // Amplitude(s) for diagram number 127
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[55], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 127 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[15] -= amp_sv[0];

    // *** DIAGRAM 128 OF 1240 ***
    // Wavefunction(s) for diagram number 128
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[6], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[57] );
    // Amplitude(s) for diagram number 128
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[57], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 128 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 129 OF 1240 ***
    // Wavefunction(s) for diagram number 129
    // (none)
    // Amplitude(s) for diagram number 129
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[57], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 129 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[21] -= amp_sv[0];

    // *** DIAGRAM 130 OF 1240 ***
    // Wavefunction(s) for diagram number 130
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[52], w_fp[34], COUPs[1], 1.0, 0., 0., w_fp[58] );
    // Amplitude(s) for diagram number 130
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[24], w_fp[6], w_fp[58], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 130 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 131 OF 1240 ***
    // Wavefunction(s) for diagram number 131
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[24], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[59] );
    // Amplitude(s) for diagram number 131
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[59], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 131 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 132 OF 1240 ***
    // Wavefunction(s) for diagram number 132
    // (none)
    // Amplitude(s) for diagram number 132
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[57], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 132 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 133 OF 1240 ***
    // Wavefunction(s) for diagram number 133
    // (none)
    // Amplitude(s) for diagram number 133
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[5], w_fp[58], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 133 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];

    // *** DIAGRAM 134 OF 1240 ***
    // Wavefunction(s) for diagram number 134
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[27], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[60] );
    // Amplitude(s) for diagram number 134
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[60], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 134 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 135 OF 1240 ***
    // Wavefunction(s) for diagram number 135
    // (none)
    // Amplitude(s) for diagram number 135
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[55], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 135 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 136 OF 1240 ***
    // Wavefunction(s) for diagram number 136
    // (none)
    // Amplitude(s) for diagram number 136
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[4], w_fp[29], w_fp[58], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 136 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 137 OF 1240 ***
    // Wavefunction(s) for diagram number 137
    // (none)
    // Amplitude(s) for diagram number 137
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[9], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 137 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 138 OF 1240 ***
    // Wavefunction(s) for diagram number 138
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[29], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[58] );
    // Amplitude(s) for diagram number 138
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[58], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 138 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 139 OF 1240 ***
    // Wavefunction(s) for diagram number 139
    // (none)
    // Amplitude(s) for diagram number 139
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[34], w_fp[30], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[34], w_fp[31], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[34], w_fp[32], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 140 OF 1240 ***
    // Wavefunction(s) for diagram number 140
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[61] );
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[34], COUPs[1], 1.0, 0., 0., w_fp[62] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[63] );
    // Amplitude(s) for diagram number 140
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[63], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 140 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 141 OF 1240 ***
    // Wavefunction(s) for diagram number 141
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[61], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[64] );
    // Amplitude(s) for diagram number 141
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[64], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 141 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 142 OF 1240 ***
    // Wavefunction(s) for diagram number 142
    // (none)
    // Amplitude(s) for diagram number 142
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 143 OF 1240 ***
    // Wavefunction(s) for diagram number 143
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[61], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[65] );
    // Amplitude(s) for diagram number 143
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[55], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 143 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 144 OF 1240 ***
    // Wavefunction(s) for diagram number 144
    // (none)
    // Amplitude(s) for diagram number 144
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[55], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 144 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];

    // *** DIAGRAM 145 OF 1240 ***
    // Wavefunction(s) for diagram number 145
    // (none)
    // Amplitude(s) for diagram number 145
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[57], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 145 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 146 OF 1240 ***
    // Wavefunction(s) for diagram number 146
    // (none)
    // Amplitude(s) for diagram number 146
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[57], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 146 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 147 OF 1240 ***
    // Wavefunction(s) for diagram number 147
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[61], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[66] );
    // Amplitude(s) for diagram number 147
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[66], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 147 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 148 OF 1240 ***
    // Wavefunction(s) for diagram number 148
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[38], w_fp[34], COUPs[1], 1.0, 0., 0., w_fp[67] );
    // Amplitude(s) for diagram number 148
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[6], w_fp[67], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 148 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];

    // *** DIAGRAM 149 OF 1240 ***
    // Wavefunction(s) for diagram number 149
    // (none)
    // Amplitude(s) for diagram number 149
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[57], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 149 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 150 OF 1240 ***
    // Wavefunction(s) for diagram number 150
    // (none)
    // Amplitude(s) for diagram number 150
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[66], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 150 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 151 OF 1240 ***
    // Wavefunction(s) for diagram number 151
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[41], w_fp[34], COUPs[1], 1.0, 0., 0., w_fp[68] );
    // Amplitude(s) for diagram number 151
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[68], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 151 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];

    // *** DIAGRAM 152 OF 1240 ***
    // Wavefunction(s) for diagram number 152
    // (none)
    // Amplitude(s) for diagram number 152
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[55], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 152 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 153 OF 1240 ***
    // Wavefunction(s) for diagram number 153
    // (none)
    // Amplitude(s) for diagram number 153
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[66], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 153 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];

    // *** DIAGRAM 154 OF 1240 ***
    // Wavefunction(s) for diagram number 154
    // (none)
    // Amplitude(s) for diagram number 154
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[29], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 154 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 155 OF 1240 ***
    // Wavefunction(s) for diagram number 155
    // (none)
    // Amplitude(s) for diagram number 155
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[58], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 155 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];

    // *** DIAGRAM 156 OF 1240 ***
    // Wavefunction(s) for diagram number 156
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[66] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[69] );
    // Amplitude(s) for diagram number 156
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[69], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 156 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 157 OF 1240 ***
    // Wavefunction(s) for diagram number 157
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[66], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[70] );
    // Amplitude(s) for diagram number 157
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[70], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 157 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 158 OF 1240 ***
    // Wavefunction(s) for diagram number 158
    // (none)
    // Amplitude(s) for diagram number 158
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 159 OF 1240 ***
    // Wavefunction(s) for diagram number 159
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[66], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[71] );
    // Amplitude(s) for diagram number 159
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[9], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 159 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 160 OF 1240 ***
    // Wavefunction(s) for diagram number 160
    // (none)
    // Amplitude(s) for diagram number 160
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[70], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 160 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];

    // *** DIAGRAM 161 OF 1240 ***
    // Wavefunction(s) for diagram number 161
    // (none)
    // Amplitude(s) for diagram number 161
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[71], w_fp[57], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 161 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 162 OF 1240 ***
    // Wavefunction(s) for diagram number 162
    // (none)
    // Amplitude(s) for diagram number 162
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[57], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 162 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];

    // *** DIAGRAM 163 OF 1240 ***
    // Wavefunction(s) for diagram number 163
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[66], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[72] );
    // Amplitude(s) for diagram number 163
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[72], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 163 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 164 OF 1240 ***
    // Wavefunction(s) for diagram number 164
    FFV1P0_3<W_ACCESS, CD_ACCESS>( w_fp[46], w_fp[34], COUPs[1], 1.0, 0., 0., w_fp[73] );
    // Amplitude(s) for diagram number 164
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[6], w_fp[73], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 164 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];

    // *** DIAGRAM 165 OF 1240 ***
    // Wavefunction(s) for diagram number 165
    // (none)
    // Amplitude(s) for diagram number 165
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[57], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 165 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 166 OF 1240 ***
    // Wavefunction(s) for diagram number 166
    // (none)
    // Amplitude(s) for diagram number 166
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[72], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 166 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 167 OF 1240 ***
    // Wavefunction(s) for diagram number 167
    // (none)
    // Amplitude(s) for diagram number 167
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[68], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 167 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];

    // *** DIAGRAM 168 OF 1240 ***
    // Wavefunction(s) for diagram number 168
    // (none)
    // Amplitude(s) for diagram number 168
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[9], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 168 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 169 OF 1240 ***
    // Wavefunction(s) for diagram number 169
    // (none)
    // Amplitude(s) for diagram number 169
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[72], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 169 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];

    // *** DIAGRAM 170 OF 1240 ***
    // Wavefunction(s) for diagram number 170
    // (none)
    // Amplitude(s) for diagram number 170
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[27], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 170 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 171 OF 1240 ***
    // Wavefunction(s) for diagram number 171
    // (none)
    // Amplitude(s) for diagram number 171
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[60], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 171 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];

    // *** DIAGRAM 172 OF 1240 ***
    // Wavefunction(s) for diagram number 172
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[1], w_fp[6], COUPs[0], 1.0, 0., 0., w_fp[72] );
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[74] );
    // Amplitude(s) for diagram number 172
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[74], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 172 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 173 OF 1240 ***
    // Wavefunction(s) for diagram number 173
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[72], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[75] );
    // Amplitude(s) for diagram number 173
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[75], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 173 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 174 OF 1240 ***
    // Wavefunction(s) for diagram number 174
    // (none)
    // Amplitude(s) for diagram number 174
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[5], w_fp[62], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 175 OF 1240 ***
    // Wavefunction(s) for diagram number 175
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[3], w_fp[72], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[76] );
    // Amplitude(s) for diagram number 175
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[76], w_fp[9], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 175 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 176 OF 1240 ***
    // Wavefunction(s) for diagram number 176
    // (none)
    // Amplitude(s) for diagram number 176
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[9], w_fp[75], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 176 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];

    // *** DIAGRAM 177 OF 1240 ***
    // Wavefunction(s) for diagram number 177
    // (none)
    // Amplitude(s) for diagram number 177
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[76], w_fp[55], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 177 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 178 OF 1240 ***
    // Wavefunction(s) for diagram number 178
    // (none)
    // Amplitude(s) for diagram number 178
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[55], w_fp[74], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 178 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];

    // *** DIAGRAM 179 OF 1240 ***
    // Wavefunction(s) for diagram number 179
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[72], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[77] );
    // Amplitude(s) for diagram number 179
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 179 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 180 OF 1240 ***
    // Wavefunction(s) for diagram number 180
    // (none)
    // Amplitude(s) for diagram number 180
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[5], w_fp[73], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 180 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];

    // *** DIAGRAM 181 OF 1240 ***
    // Wavefunction(s) for diagram number 181
    // (none)
    // Amplitude(s) for diagram number 181
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[55], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 181 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 182 OF 1240 ***
    // Wavefunction(s) for diagram number 182
    // (none)
    // Amplitude(s) for diagram number 182
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[77], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 182 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 183 OF 1240 ***
    // Wavefunction(s) for diagram number 183
    // (none)
    // Amplitude(s) for diagram number 183
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[4], w_fp[67], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 183 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];

    // *** DIAGRAM 184 OF 1240 ***
    // Wavefunction(s) for diagram number 184
    // (none)
    // Amplitude(s) for diagram number 184
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[9], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 184 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 185 OF 1240 ***
    // Wavefunction(s) for diagram number 185
    // (none)
    // Amplitude(s) for diagram number 185
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[77], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 185 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];

    // *** DIAGRAM 186 OF 1240 ***
    // Wavefunction(s) for diagram number 186
    // (none)
    // Amplitude(s) for diagram number 186
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[72], w_fp[24], w_fp[62], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 186 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 187 OF 1240 ***
    // Wavefunction(s) for diagram number 187
    // (none)
    // Amplitude(s) for diagram number 187
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[59], w_fp[72], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 187 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];

    // *** DIAGRAM 188 OF 1240 ***
    // Wavefunction(s) for diagram number 188
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[34], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[77] );
    // Amplitude(s) for diagram number 188
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 188 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[5] -= amp_sv[0];

    // *** DIAGRAM 189 OF 1240 ***
    // Wavefunction(s) for diagram number 189
    // (none)
    // Amplitude(s) for diagram number 189
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[77], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 189 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] -= amp_sv[0];

    // *** DIAGRAM 190 OF 1240 ***
    // Wavefunction(s) for diagram number 190
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[46], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[78] );
    // Amplitude(s) for diagram number 190
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[78], w_fp[55], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 190 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] -= amp_sv[0];

    // *** DIAGRAM 191 OF 1240 ***
    // Wavefunction(s) for diagram number 191
    // (none)
    // Amplitude(s) for diagram number 191
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[53], w_fp[55], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 191 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[13] -= amp_sv[0];

    // *** DIAGRAM 192 OF 1240 ***
    // Wavefunction(s) for diagram number 192
    // (none)
    // Amplitude(s) for diagram number 192
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[78], w_fp[57], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 192 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[22] -= amp_sv[0];

    // *** DIAGRAM 193 OF 1240 ***
    // Wavefunction(s) for diagram number 193
    // (none)
    // Amplitude(s) for diagram number 193
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[7], w_fp[57], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 193 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[19] -= amp_sv[0];

    // *** DIAGRAM 194 OF 1240 ***
    // Wavefunction(s) for diagram number 194
    // (none)
    // Amplitude(s) for diagram number 194
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[77], w_fp[29], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 194 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 195 OF 1240 ***
    // Wavefunction(s) for diagram number 195
    // (none)
    // Amplitude(s) for diagram number 195
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[29], w_fp[73], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 195 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];

    // *** DIAGRAM 196 OF 1240 ***
    // Wavefunction(s) for diagram number 196
    // (none)
    // Amplitude(s) for diagram number 196
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[46], w_fp[58], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 196 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 197 OF 1240 ***
    // Wavefunction(s) for diagram number 197
    // (none)
    // Amplitude(s) for diagram number 197
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[25], w_fp[77], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 197 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[4] -= amp_sv[0];

    // *** DIAGRAM 198 OF 1240 ***
    // Wavefunction(s) for diagram number 198
    // (none)
    // Amplitude(s) for diagram number 198
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[77], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 198 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= amp_sv[0];

    // *** DIAGRAM 199 OF 1240 ***
    // Wavefunction(s) for diagram number 199
    FFV1_2<W_ACCESS, CD_ACCESS>( w_fp[38], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[58] );
    // Amplitude(s) for diagram number 199
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[58], w_fp[9], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 199 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[10] -= amp_sv[0];

    // *** DIAGRAM 200 OF 1240 ***
    // Wavefunction(s) for diagram number 200
    // (none)
    // Amplitude(s) for diagram number 200
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[9], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 200 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[7] -= amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcol( jamps, icol ) += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 9 );
    storeWf( wfs, w_cx, nevt, 12 );
    storeWf( wfs, w_cx, nevt, 14 );
    storeWf( wfs, w_cx, nevt, 15 );
    storeWf( wfs, w_cx, nevt, 16 );
    storeWf( wfs, w_cx, nevt, 17 );
    storeWf( wfs, w_cx, nevt, 18 );
    storeWf( wfs, w_cx, nevt, 19 );
    storeWf( wfs, w_cx, nevt, 21 );
    storeWf( wfs, w_cx, nevt, 22 );
    storeWf( wfs, w_cx, nevt, 26 );
    storeWf( wfs, w_cx, nevt, 34 );
    storeWf( wfs, w_cx, nevt, 42 );
    storeWf( wfs, w_cx, nevt, 52 );
    storeWf( wfs, w_cx, nevt, 55 );
    storeWf( wfs, w_cx, nevt, 56 );
    storeWf( wfs, w_cx, nevt, 57 );
    storeWf( wfs, w_cx, nevt, 58 );
    storeWf( wfs, w_cx, nevt, 59 );
    storeWf( wfs, w_cx, nevt, 60 );
    storeWf( wfs, w_cx, nevt, 61 );
    storeWf( wfs, w_cx, nevt, 62 );
    storeWf( wfs, w_cx, nevt, 63 );
    storeWf( wfs, w_cx, nevt, 64 );
    storeWf( wfs, w_cx, nevt, 65 );
    storeWf( wfs, w_cx, nevt, 66 );
    storeWf( wfs, w_cx, nevt, 67 );
    storeWf( wfs, w_cx, nevt, 68 );
    storeWf( wfs, w_cx, nevt, 69 );
    storeWf( wfs, w_cx, nevt, 70 );
    storeWf( wfs, w_cx, nevt, 71 );
    storeWf( wfs, w_cx, nevt, 72 );
    storeWf( wfs, w_cx, nevt, 73 );
    storeWf( wfs, w_cx, nevt, 74 );
    storeWf( wfs, w_cx, nevt, 75 );
    storeWf( wfs, w_cx, nevt, 76 );
    storeWf( wfs, w_cx, nevt, 77 );
    storeWf( wfs, w_cx, nevt, 78 );
#endif
  }

  //--------------------------------------------------------------------------
}
