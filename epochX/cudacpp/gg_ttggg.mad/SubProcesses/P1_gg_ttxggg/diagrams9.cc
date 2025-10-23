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
  diagramgroup9( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 22 );
    retrieveWf( wfs, w_cx, nevt, 23 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 31 );
    retrieveWf( wfs, w_cx, nevt, 32 );
    retrieveWf( wfs, w_cx, nevt, 34 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 41 );
    retrieveWf( wfs, w_cx, nevt, 42 );
    retrieveWf( wfs, w_cx, nevt, 43 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 48 );
    retrieveWf( wfs, w_cx, nevt, 49 );
    retrieveWf( wfs, w_cx, nevt, 50 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 93 );
    retrieveWf( wfs, w_cx, nevt, 96 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 106 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 112 );
#endif

    // *** DIAGRAM 801 OF 1240 ***
    // Wavefunction(s) for diagram number 801
    // (none)
    // Amplitude(s) for diagram number 801
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[102], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 801 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[75] -= amp_sv[0];

    // *** DIAGRAM 802 OF 1240 ***
    // Wavefunction(s) for diagram number 802
    // (none)
    // Amplitude(s) for diagram number 802
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[44], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 802 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[91] -= amp_sv[0];

    // *** DIAGRAM 803 OF 1240 ***
    // Wavefunction(s) for diagram number 803
    // (none)
    // Amplitude(s) for diagram number 803
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[43], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 803 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[85] -= amp_sv[0];

    // *** DIAGRAM 804 OF 1240 ***
    // Wavefunction(s) for diagram number 804
    // (none)
    // Amplitude(s) for diagram number 804
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[39], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 804 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 805 OF 1240 ***
    // Wavefunction(s) for diagram number 805
    // (none)
    // Amplitude(s) for diagram number 805
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[102], w_fp[27], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 805 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 806 OF 1240 ***
    // Wavefunction(s) for diagram number 806
    // (none)
    // Amplitude(s) for diagram number 806
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[39], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 806 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];

    // *** DIAGRAM 807 OF 1240 ***
    // Wavefunction(s) for diagram number 807
    // (none)
    // Amplitude(s) for diagram number 807
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[49], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 807 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[112] -= amp_sv[0];

    // *** DIAGRAM 808 OF 1240 ***
    // Wavefunction(s) for diagram number 808
    // (none)
    // Amplitude(s) for diagram number 808
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[50], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 808 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 809 OF 1240 ***
    // Wavefunction(s) for diagram number 809
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[47], w_fp[0], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[113] );
    // Amplitude(s) for diagram number 809
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[113], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 809 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 810 OF 1240 ***
    // Wavefunction(s) for diagram number 810
    // (none)
    // Amplitude(s) for diagram number 810
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[113], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 810 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[99] -= amp_sv[0];

    // *** DIAGRAM 811 OF 1240 ***
    // Wavefunction(s) for diagram number 811
    // (none)
    // Amplitude(s) for diagram number 811
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[50], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 811 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[115] -= amp_sv[0];

    // *** DIAGRAM 812 OF 1240 ***
    // Wavefunction(s) for diagram number 812
    // (none)
    // Amplitude(s) for diagram number 812
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[49], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 812 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[109] -= amp_sv[0];

    // *** DIAGRAM 813 OF 1240 ***
    // Wavefunction(s) for diagram number 813
    // (none)
    // Amplitude(s) for diagram number 813
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[47], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 813 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 814 OF 1240 ***
    // Wavefunction(s) for diagram number 814
    // (none)
    // Amplitude(s) for diagram number 814
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[113], w_fp[24], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 814 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 815 OF 1240 ***
    // Wavefunction(s) for diagram number 815
    // (none)
    // Amplitude(s) for diagram number 815
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[47], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 815 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 816 OF 1240 ***
    // Wavefunction(s) for diagram number 816
    // (none)
    // Amplitude(s) for diagram number 816
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[17], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 816 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 817 OF 1240 ***
    // Wavefunction(s) for diagram number 817
    // (none)
    // Amplitude(s) for diagram number 817
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[2], w_fp[42], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 817 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 818 OF 1240 ***
    // Wavefunction(s) for diagram number 818
    // (none)
    // Amplitude(s) for diagram number 818
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[98], w_fp[103], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 818 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 819 OF 1240 ***
    // Wavefunction(s) for diagram number 819
    // (none)
    // Amplitude(s) for diagram number 819
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[2], w_fp[98], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 819 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];

    // *** DIAGRAM 820 OF 1240 ***
    // Wavefunction(s) for diagram number 820
    // (none)
    // Amplitude(s) for diagram number 820
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[103], w_fp[42], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 820 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 821 OF 1240 ***
    // Wavefunction(s) for diagram number 821
    // (none)
    // Amplitude(s) for diagram number 821
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[21], w_fp[17], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 821 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 822 OF 1240 ***
    // Wavefunction(s) for diagram number 822
    // (none)
    // Amplitude(s) for diagram number 822
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[85], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[112], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[111], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 823 OF 1240 ***
    // Wavefunction(s) for diagram number 823
    // (none)
    // Amplitude(s) for diagram number 823
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[15], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 823 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 824 OF 1240 ***
    // Wavefunction(s) for diagram number 824
    // (none)
    // Amplitude(s) for diagram number 824
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[2], w_fp[16], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 824 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];

    // *** DIAGRAM 825 OF 1240 ***
    // Wavefunction(s) for diagram number 825
    // (none)
    // Amplitude(s) for diagram number 825
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[101], w_fp[103], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 825 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 826 OF 1240 ***
    // Wavefunction(s) for diagram number 826
    // (none)
    // Amplitude(s) for diagram number 826
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[2], w_fp[101], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 826 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 827 OF 1240 ***
    // Wavefunction(s) for diagram number 827
    // (none)
    // Amplitude(s) for diagram number 827
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[103], w_fp[16], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 827 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 828 OF 1240 ***
    // Wavefunction(s) for diagram number 828
    // (none)
    // Amplitude(s) for diagram number 828
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[22], w_fp[15], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 828 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 829 OF 1240 ***
    // Wavefunction(s) for diagram number 829
    // (none)
    // Amplitude(s) for diagram number 829
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[9], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[110], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[109], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 830 OF 1240 ***
    // Wavefunction(s) for diagram number 830
    // (none)
    // Amplitude(s) for diagram number 830
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[18], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 830 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 831 OF 1240 ***
    // Wavefunction(s) for diagram number 831
    // (none)
    // Amplitude(s) for diagram number 831
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[2], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 831 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 832 OF 1240 ***
    // Wavefunction(s) for diagram number 832
    // (none)
    // Amplitude(s) for diagram number 832
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[96], w_fp[103], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 832 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 833 OF 1240 ***
    // Wavefunction(s) for diagram number 833
    // (none)
    // Amplitude(s) for diagram number 833
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[2], w_fp[96], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 833 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];

    // *** DIAGRAM 834 OF 1240 ***
    // Wavefunction(s) for diagram number 834
    // (none)
    // Amplitude(s) for diagram number 834
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[103], w_fp[19], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 834 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 835 OF 1240 ***
    // Wavefunction(s) for diagram number 835
    // (none)
    // Amplitude(s) for diagram number 835
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[18], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 835 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 836 OF 1240 ***
    // Wavefunction(s) for diagram number 836
    // (none)
    // Amplitude(s) for diagram number 836
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[34], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 837 OF 1240 ***
    // Wavefunction(s) for diagram number 837
    // (none)
    // Amplitude(s) for diagram number 837
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[2], w_fp[30], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[2], w_fp[31], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[2], w_fp[32], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 838 OF 1240 ***
    // Wavefunction(s) for diagram number 838
    // (none)
    // Amplitude(s) for diagram number 838
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[92], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[52], w_fp[2], w_fp[106], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 839 OF 1240 ***
    // Wavefunction(s) for diagram number 839
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], COUPs[0], 1.0, 0., 0., w_fp[90] );
    // Amplitude(s) for diagram number 839
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[10], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 839 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 840 OF 1240 ***
    // Wavefunction(s) for diagram number 840
    // (none)
    // Amplitude(s) for diagram number 840
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[11], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 840 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 841 OF 1240 ***
    // Wavefunction(s) for diagram number 841
    // (none)
    // Amplitude(s) for diagram number 841
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[6], w_fp[90], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[6], w_fp[90], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[6], w_fp[90], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 842 OF 1240 ***
    // Wavefunction(s) for diagram number 842
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], COUPs[0], 1.0, 0., 0., w_fp[56] );
    // Amplitude(s) for diagram number 842
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[63], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 842 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 843 OF 1240 ***
    // Wavefunction(s) for diagram number 843
    // (none)
    // Amplitude(s) for diagram number 843
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[64], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 843 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];

    // *** DIAGRAM 844 OF 1240 ***
    // Wavefunction(s) for diagram number 844
    // (none)
    // Amplitude(s) for diagram number 844
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[6], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 845 OF 1240 ***
    // Wavefunction(s) for diagram number 845
    // (none)
    // Amplitude(s) for diagram number 845
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[63], w_fp[11], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 845 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 846 OF 1240 ***
    // Wavefunction(s) for diagram number 846
    // (none)
    // Amplitude(s) for diagram number 846
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[64], w_fp[10], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 846 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];

    // *** DIAGRAM 847 OF 1240 ***
    // Wavefunction(s) for diagram number 847
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[103] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[22] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[21] );
    // Amplitude(s) for diagram number 847
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[103], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[22], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[6], w_fp[21], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 848 OF 1240 ***
    // Wavefunction(s) for diagram number 848
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[105] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[95] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[107] );
    // Amplitude(s) for diagram number 848
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[105], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[95], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[5], w_fp[107], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 849 OF 1240 ***
    // Wavefunction(s) for diagram number 849
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[115] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[116] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[5], COUPs[2], 1.0, 0., 0., w_fp[117] );
    // Amplitude(s) for diagram number 849
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[6], w_fp[115], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[6], w_fp[116], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[6], w_fp[117], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];

    // *** DIAGRAM 850 OF 1240 ***
    // Wavefunction(s) for diagram number 850
    VVVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[118] );
    VVVV3P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[119] );
    VVVV4P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[8], w_fp[6], COUPs[2], 1.0, 0., 0., w_fp[120] );
    // Amplitude(s) for diagram number 850
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[118], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[119], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[5], w_fp[120], COUPs[0], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 851 OF 1240 ***
    // Wavefunction(s) for diagram number 851
    // (none)
    // Amplitude(s) for diagram number 851
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[8], w_fp[29], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[8], w_fp[29], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[61], w_fp[8], w_fp[29], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];

    // *** DIAGRAM 852 OF 1240 ***
    // Wavefunction(s) for diagram number 852
    // (none)
    // Amplitude(s) for diagram number 852
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[29], w_fp[90], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 852 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
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

    // *** DIAGRAM 853 OF 1240 ***
    // Wavefunction(s) for diagram number 853
    // (none)
    // Amplitude(s) for diagram number 853
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[29], w_fp[56], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 853 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 854 OF 1240 ***
    // Wavefunction(s) for diagram number 854
    // (none)
    // Amplitude(s) for diagram number 854
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[61], w_fp[8], w_fp[96], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 854 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];

    // *** DIAGRAM 855 OF 1240 ***
    // Wavefunction(s) for diagram number 855
    // (none)
    // Amplitude(s) for diagram number 855
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[45], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 855 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 856 OF 1240 ***
    // Wavefunction(s) for diagram number 856
    // (none)
    // Amplitude(s) for diagram number 856
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[44], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 856 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];

    // *** DIAGRAM 857 OF 1240 ***
    // Wavefunction(s) for diagram number 857
    // (none)
    // Amplitude(s) for diagram number 857
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[102], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 857 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 858 OF 1240 ***
    // Wavefunction(s) for diagram number 858
    // (none)
    // Amplitude(s) for diagram number 858
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[102], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 858 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];

    // *** DIAGRAM 859 OF 1240 ***
    // Wavefunction(s) for diagram number 859
    // (none)
    // Amplitude(s) for diagram number 859
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[44], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 859 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 860 OF 1240 ***
    // Wavefunction(s) for diagram number 860
    // (none)
    // Amplitude(s) for diagram number 860
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[64], w_fp[45], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 860 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 861 OF 1240 ***
    // Wavefunction(s) for diagram number 861
    // (none)
    // Amplitude(s) for diagram number 861
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[105], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[95], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[39], w_fp[107], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 862 OF 1240 ***
    // Wavefunction(s) for diagram number 862
    // (none)
    // Amplitude(s) for diagram number 862
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[39], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 862 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 863 OF 1240 ***
    // Wavefunction(s) for diagram number 863
    // (none)
    // Amplitude(s) for diagram number 863
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[102], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 863 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 864 OF 1240 ***
    // Wavefunction(s) for diagram number 864
    // (none)
    // Amplitude(s) for diagram number 864
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[39], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 864 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 865 OF 1240 ***
    // Wavefunction(s) for diagram number 865
    // (none)
    // Amplitude(s) for diagram number 865
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[51], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 865 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 866 OF 1240 ***
    // Wavefunction(s) for diagram number 866
    // (none)
    // Amplitude(s) for diagram number 866
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[50], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 866 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 867 OF 1240 ***
    // Wavefunction(s) for diagram number 867
    // (none)
    // Amplitude(s) for diagram number 867
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[113], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 867 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 868 OF 1240 ***
    // Wavefunction(s) for diagram number 868
    // (none)
    // Amplitude(s) for diagram number 868
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[113], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 868 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];

    // *** DIAGRAM 869 OF 1240 ***
    // Wavefunction(s) for diagram number 869
    // (none)
    // Amplitude(s) for diagram number 869
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[50], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 869 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 870 OF 1240 ***
    // Wavefunction(s) for diagram number 870
    // (none)
    // Amplitude(s) for diagram number 870
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[63], w_fp[51], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 870 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 871 OF 1240 ***
    // Wavefunction(s) for diagram number 871
    // (none)
    // Amplitude(s) for diagram number 871
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 872 OF 1240 ***
    // Wavefunction(s) for diagram number 872
    // (none)
    // Amplitude(s) for diagram number 872
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[47], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 872 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 873 OF 1240 ***
    // Wavefunction(s) for diagram number 873
    // (none)
    // Amplitude(s) for diagram number 873
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[113], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 873 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 874 OF 1240 ***
    // Wavefunction(s) for diagram number 874
    // (none)
    // Amplitude(s) for diagram number 874
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[47], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 874 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 875 OF 1240 ***
    // Wavefunction(s) for diagram number 875
    // (none)
    // Amplitude(s) for diagram number 875
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[23], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 875 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 876 OF 1240 ***
    // Wavefunction(s) for diagram number 876
    // (none)
    // Amplitude(s) for diagram number 876
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[2], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 876 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];

    // *** DIAGRAM 877 OF 1240 ***
    // Wavefunction(s) for diagram number 877
    // (none)
    // Amplitude(s) for diagram number 877
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[93], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 877 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 878 OF 1240 ***
    // Wavefunction(s) for diagram number 878
    // (none)
    // Amplitude(s) for diagram number 878
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[104], w_fp[2], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 878 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 879 OF 1240 ***
    // Wavefunction(s) for diagram number 879
    // (none)
    // Amplitude(s) for diagram number 879
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[48], w_fp[93], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 879 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 880 OF 1240 ***
    // Wavefunction(s) for diagram number 880
    // (none)
    // Amplitude(s) for diagram number 880
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[64], w_fp[23], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 880 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 881 OF 1240 ***
    // Wavefunction(s) for diagram number 881
    // (none)
    // Amplitude(s) for diagram number 881
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[105], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[95], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[38], w_fp[2], w_fp[107], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 882 OF 1240 ***
    // Wavefunction(s) for diagram number 882
    // (none)
    // Amplitude(s) for diagram number 882
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[20], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 882 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 883 OF 1240 ***
    // Wavefunction(s) for diagram number 883
    // (none)
    // Amplitude(s) for diagram number 883
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[2], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 883 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];

    // *** DIAGRAM 884 OF 1240 ***
    // Wavefunction(s) for diagram number 884
    // (none)
    // Amplitude(s) for diagram number 884
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[93], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 884 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 885 OF 1240 ***
    // Wavefunction(s) for diagram number 885
    // (none)
    // Amplitude(s) for diagram number 885
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[2], w_fp[63], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 885 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 886 OF 1240 ***
    // Wavefunction(s) for diagram number 886
    // (none)
    // Amplitude(s) for diagram number 886
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[93], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 886 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 887 OF 1240 ***
    // Wavefunction(s) for diagram number 887
    // (none)
    // Amplitude(s) for diagram number 887
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[63], w_fp[20], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 887 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 888 OF 1240 ***
    // Wavefunction(s) for diagram number 888
    // (none)
    // Amplitude(s) for diagram number 888
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[22], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[41], w_fp[2], w_fp[21], COUPs[1], 1.0, &amp_fp[0] );
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 889 OF 1240 ***
    // Wavefunction(s) for diagram number 889
    // (none)
    // Amplitude(s) for diagram number 889
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[18], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 889 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 890 OF 1240 ***
    // Wavefunction(s) for diagram number 890
    // (none)
    // Amplitude(s) for diagram number 890
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 890 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 891 OF 1240 ***
    // Wavefunction(s) for diagram number 891
    // (none)
    // Amplitude(s) for diagram number 891
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[93], w_fp[96], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 891 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 892 OF 1240 ***
    // Wavefunction(s) for diagram number 892
    // (none)
    // Amplitude(s) for diagram number 892
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[2], w_fp[96], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 892 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 893 OF 1240 ***
    // Wavefunction(s) for diagram number 893
    // (none)
    // Amplitude(s) for diagram number 893
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[93], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 893 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];

    // *** DIAGRAM 894 OF 1240 ***
    // Wavefunction(s) for diagram number 894
    // (none)
    // Amplitude(s) for diagram number 894
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[18], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 894 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];

    // *** DIAGRAM 895 OF 1240 ***
    // Wavefunction(s) for diagram number 895
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[66], COUPs[0], 1.0, 0., 0., w_fp[65] );
    // Amplitude(s) for diagram number 895
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[13], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 895 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
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

    // *** DIAGRAM 896 OF 1240 ***
    // Wavefunction(s) for diagram number 896
    // (none)
    // Amplitude(s) for diagram number 896
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[65], w_fp[11], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 896 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
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

    // *** DIAGRAM 897 OF 1240 ***
    // Wavefunction(s) for diagram number 897
    // (none)
    // Amplitude(s) for diagram number 897
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[6], w_fp[65], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[6], w_fp[65], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[8], w_fp[4], w_fp[6], w_fp[65], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 898 OF 1240 ***
    // Wavefunction(s) for diagram number 898
    // (none)
    // Amplitude(s) for diagram number 898
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[69], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 898 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
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

    // *** DIAGRAM 899 OF 1240 ***
    // Wavefunction(s) for diagram number 899
    // (none)
    // Amplitude(s) for diagram number 899
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[56], w_fp[70], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    if( channelId == 899 ) numerators_sv += cxabs2( amp_sv[0] );
    if( channelId != 0 ) denominators_sv += cxabs2( amp_sv[0] );
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

    // *** DIAGRAM 900 OF 1240 ***
    // Wavefunction(s) for diagram number 900
    // (none)
    // Amplitude(s) for diagram number 900
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[66], w_fp[4], w_fp[6], w_fp[56], COUPs[2], 1.0, &amp_fp[0] );
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 56 );
    storeWf( wfs, w_cx, nevt, 65 );
    storeWf( wfs, w_cx, nevt, 90 );
    storeWf( wfs, w_cx, nevt, 95 );
    storeWf( wfs, w_cx, nevt, 103 );
    storeWf( wfs, w_cx, nevt, 105 );
    storeWf( wfs, w_cx, nevt, 107 );
    storeWf( wfs, w_cx, nevt, 113 );
    storeWf( wfs, w_cx, nevt, 115 );
    storeWf( wfs, w_cx, nevt, 116 );
    storeWf( wfs, w_cx, nevt, 117 );
    storeWf( wfs, w_cx, nevt, 118 );
    storeWf( wfs, w_cx, nevt, 119 );
    storeWf( wfs, w_cx, nevt, 120 );
#endif
  }

  //--------------------------------------------------------------------------
}
