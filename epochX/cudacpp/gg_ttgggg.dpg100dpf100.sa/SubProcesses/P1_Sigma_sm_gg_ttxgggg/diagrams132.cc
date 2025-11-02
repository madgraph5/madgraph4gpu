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

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup132( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                   fptype* jamps,                  // output jamps[ncolor*2*nevt]
                   const int nGoodHel,             // input: number of good helicities
                   const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                   cxtype_sv* jamps,               // output jamps[ncolor*2*neppV]
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
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** RETRIEVE WAVEFUNCTIONS FROM PREVIOUS DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) retrieveWf( wfs, w_cx, nevt, iwf );
    retrieveWf( wfs, w_cx, nevt, 0 );
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 78 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 107 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 126 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 192 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 253 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 261 );
    retrieveWf( wfs, w_cx, nevt, 278 );
    retrieveWf( wfs, w_cx, nevt, 280 );
    retrieveWf( wfs, w_cx, nevt, 352 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 364 );
    retrieveWf( wfs, w_cx, nevt, 365 );
    retrieveWf( wfs, w_cx, nevt, 370 );
    retrieveWf( wfs, w_cx, nevt, 431 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 489 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 598 );
    retrieveWf( wfs, w_cx, nevt, 599 );
    retrieveWf( wfs, w_cx, nevt, 602 );
    retrieveWf( wfs, w_cx, nevt, 603 );
    retrieveWf( wfs, w_cx, nevt, 611 );
    retrieveWf( wfs, w_cx, nevt, 612 );
    retrieveWf( wfs, w_cx, nevt, 614 );
    retrieveWf( wfs, w_cx, nevt, 615 );
    retrieveWf( wfs, w_cx, nevt, 617 );
    retrieveWf( wfs, w_cx, nevt, 624 );
    retrieveWf( wfs, w_cx, nevt, 625 );
    retrieveWf( wfs, w_cx, nevt, 627 );
    retrieveWf( wfs, w_cx, nevt, 673 );
    retrieveWf( wfs, w_cx, nevt, 684 );
    retrieveWf( wfs, w_cx, nevt, 685 );
    retrieveWf( wfs, w_cx, nevt, 686 );
    retrieveWf( wfs, w_cx, nevt, 687 );
    retrieveWf( wfs, w_cx, nevt, 688 );
    retrieveWf( wfs, w_cx, nevt, 689 );
    retrieveWf( wfs, w_cx, nevt, 745 );
#endif
#endif

    // *** DIAGRAM 13101 OF 15495 ***
    // Wavefunction(s) for diagram number 13101
    // (none)
    // Amplitude(s) for diagram number 13101
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[12], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13102 OF 15495 ***
    // Wavefunction(s) for diagram number 13102
    // (none)
    // Amplitude(s) for diagram number 13102
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[685], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13103 OF 15495 ***
    // Wavefunction(s) for diagram number 13103
    // (none)
    // Amplitude(s) for diagram number 13103
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];

    // *** DIAGRAM 13104 OF 15495 ***
    // Wavefunction(s) for diagram number 13104
    // (none)
    // Amplitude(s) for diagram number 13104
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[69], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13105 OF 15495 ***
    // Wavefunction(s) for diagram number 13105
    // (none)
    // Amplitude(s) for diagram number 13105
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[253], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[364], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[365], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13106 OF 15495 ***
    // Wavefunction(s) for diagram number 13106
    // (none)
    // Amplitude(s) for diagram number 13106
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[599], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];

    // *** DIAGRAM 13107 OF 15495 ***
    // Wavefunction(s) for diagram number 13107
    // (none)
    // Amplitude(s) for diagram number 13107
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[598], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13108 OF 15495 ***
    // Wavefunction(s) for diagram number 13108
    // (none)
    // Amplitude(s) for diagram number 13108
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 13109 OF 15495 ***
    // Wavefunction(s) for diagram number 13109
    // (none)
    // Amplitude(s) for diagram number 13109
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[7], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 13110 OF 15495 ***
    // Wavefunction(s) for diagram number 13110
    // (none)
    // Amplitude(s) for diagram number 13110
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[214], w_fp[598], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 13111 OF 15495 ***
    // Wavefunction(s) for diagram number 13111
    // (none)
    // Amplitude(s) for diagram number 13111
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13112 OF 15495 ***
    // Wavefunction(s) for diagram number 13112
    // (none)
    // Amplitude(s) for diagram number 13112
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[599], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];

    // *** DIAGRAM 13113 OF 15495 ***
    // Wavefunction(s) for diagram number 13113
    // (none)
    // Amplitude(s) for diagram number 13113
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[40], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13114 OF 15495 ***
    // Wavefunction(s) for diagram number 13114
    // (none)
    // Amplitude(s) for diagram number 13114
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[603], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13115 OF 15495 ***
    // Wavefunction(s) for diagram number 13115
    // (none)
    // Amplitude(s) for diagram number 13115
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[40], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 13116 OF 15495 ***
    // Wavefunction(s) for diagram number 13116
    // (none)
    // Amplitude(s) for diagram number 13116
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[602], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13117 OF 15495 ***
    // Wavefunction(s) for diagram number 13117
    // (none)
    // Amplitude(s) for diagram number 13117
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[214], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 13118 OF 15495 ***
    // Wavefunction(s) for diagram number 13118
    // (none)
    // Amplitude(s) for diagram number 13118
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[7], w_fp[182], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 13119 OF 15495 ***
    // Wavefunction(s) for diagram number 13119
    // (none)
    // Amplitude(s) for diagram number 13119
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[254], w_fp[7], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 13120 OF 15495 ***
    // Wavefunction(s) for diagram number 13120
    // (none)
    // Amplitude(s) for diagram number 13120
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[182], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13121 OF 15495 ***
    // Wavefunction(s) for diagram number 13121
    // (none)
    // Amplitude(s) for diagram number 13121
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[688], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 13122 OF 15495 ***
    // Wavefunction(s) for diagram number 13122
    // (none)
    // Amplitude(s) for diagram number 13122
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[69], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[69], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[69], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 13123 OF 15495 ***
    // Wavefunction(s) for diagram number 13123
    // (none)
    // Amplitude(s) for diagram number 13123
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[69], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 13124 OF 15495 ***
    // Wavefunction(s) for diagram number 13124
    // (none)
    // Amplitude(s) for diagram number 13124
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[214], w_fp[602], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 13125 OF 15495 ***
    // Wavefunction(s) for diagram number 13125
    // (none)
    // Amplitude(s) for diagram number 13125
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[688], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13126 OF 15495 ***
    // Wavefunction(s) for diagram number 13126
    // (none)
    // Amplitude(s) for diagram number 13126
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[603], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13127 OF 15495 ***
    // Wavefunction(s) for diagram number 13127
    // (none)
    // Amplitude(s) for diagram number 13127
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[136], w_fp[214], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[137], w_fp[214], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[551], w_fp[214], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 13128 OF 15495 ***
    // Wavefunction(s) for diagram number 13128
    // (none)
    // Amplitude(s) for diagram number 13128
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[136], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[137], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[213], w_fp[551], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13129 OF 15495 ***
    // Wavefunction(s) for diagram number 13129
    // (none)
    // Amplitude(s) for diagram number 13129
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[526], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[524], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13130 OF 15495 ***
    // Wavefunction(s) for diagram number 13130
    // (none)
    // Amplitude(s) for diagram number 13130
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[521], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[526], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[524], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 13131 OF 15495 ***
    // Wavefunction(s) for diagram number 13131
    // (none)
    // Amplitude(s) for diagram number 13131
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[253], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[364], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[365], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 13132 OF 15495 ***
    // Wavefunction(s) for diagram number 13132
    // (none)
    // Amplitude(s) for diagram number 13132
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[363], w_fp[5], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13133 OF 15495 ***
    // Wavefunction(s) for diagram number 13133
    // (none)
    // Amplitude(s) for diagram number 13133
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[684], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];

    // *** DIAGRAM 13134 OF 15495 ***
    // Wavefunction(s) for diagram number 13134
    // (none)
    // Amplitude(s) for diagram number 13134
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[12], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13135 OF 15495 ***
    // Wavefunction(s) for diagram number 13135
    // (none)
    // Amplitude(s) for diagram number 13135
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[684], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13136 OF 15495 ***
    // Wavefunction(s) for diagram number 13136
    // (none)
    // Amplitude(s) for diagram number 13136
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];

    // *** DIAGRAM 13137 OF 15495 ***
    // Wavefunction(s) for diagram number 13137
    // (none)
    // Amplitude(s) for diagram number 13137
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[103], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13138 OF 15495 ***
    // Wavefunction(s) for diagram number 13138
    // (none)
    // Amplitude(s) for diagram number 13138
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[107], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[78], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[370], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13139 OF 15495 ***
    // Wavefunction(s) for diagram number 13139
    // (none)
    // Amplitude(s) for diagram number 13139
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[614], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[516] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];

    // *** DIAGRAM 13140 OF 15495 ***
    // Wavefunction(s) for diagram number 13140
    // (none)
    // Amplitude(s) for diagram number 13140
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[612], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13141 OF 15495 ***
    // Wavefunction(s) for diagram number 13141
    // (none)
    // Amplitude(s) for diagram number 13141
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];

    // *** DIAGRAM 13142 OF 15495 ***
    // Wavefunction(s) for diagram number 13142
    // (none)
    // Amplitude(s) for diagram number 13142
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[5], w_fp[611], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];

    // *** DIAGRAM 13143 OF 15495 ***
    // Wavefunction(s) for diagram number 13143
    // (none)
    // Amplitude(s) for diagram number 13143
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[214], w_fp[612], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];

    // *** DIAGRAM 13144 OF 15495 ***
    // Wavefunction(s) for diagram number 13144
    // (none)
    // Amplitude(s) for diagram number 13144
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[611], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13145 OF 15495 ***
    // Wavefunction(s) for diagram number 13145
    // (none)
    // Amplitude(s) for diagram number 13145
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[614], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];

    // *** DIAGRAM 13146 OF 15495 ***
    // Wavefunction(s) for diagram number 13146
    // (none)
    // Amplitude(s) for diagram number 13146
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[40], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13147 OF 15495 ***
    // Wavefunction(s) for diagram number 13147
    // (none)
    // Amplitude(s) for diagram number 13147
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[617], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13148 OF 15495 ***
    // Wavefunction(s) for diagram number 13148
    // (none)
    // Amplitude(s) for diagram number 13148
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[40], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];

    // *** DIAGRAM 13149 OF 15495 ***
    // Wavefunction(s) for diagram number 13149
    // (none)
    // Amplitude(s) for diagram number 13149
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[615], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13150 OF 15495 ***
    // Wavefunction(s) for diagram number 13150
    // (none)
    // Amplitude(s) for diagram number 13150
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[214], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 13151 OF 15495 ***
    // Wavefunction(s) for diagram number 13151
    // (none)
    // Amplitude(s) for diagram number 13151
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[5], w_fp[61], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 13152 OF 15495 ***
    // Wavefunction(s) for diagram number 13152
    // (none)
    // Amplitude(s) for diagram number 13152
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[363], w_fp[5], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 13153 OF 15495 ***
    // Wavefunction(s) for diagram number 13153
    // (none)
    // Amplitude(s) for diagram number 13153
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[61], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13154 OF 15495 ***
    // Wavefunction(s) for diagram number 13154
    // (none)
    // Amplitude(s) for diagram number 13154
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[687], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];

    // *** DIAGRAM 13155 OF 15495 ***
    // Wavefunction(s) for diagram number 13155
    // (none)
    // Amplitude(s) for diagram number 13155
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[103], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[103], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[103], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];

    // *** DIAGRAM 13156 OF 15495 ***
    // Wavefunction(s) for diagram number 13156
    // (none)
    // Amplitude(s) for diagram number 13156
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[103], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];

    // *** DIAGRAM 13157 OF 15495 ***
    // Wavefunction(s) for diagram number 13157
    // (none)
    // Amplitude(s) for diagram number 13157
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[214], w_fp[615], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];

    // *** DIAGRAM 13158 OF 15495 ***
    // Wavefunction(s) for diagram number 13158
    // (none)
    // Amplitude(s) for diagram number 13158
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[687], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13159 OF 15495 ***
    // Wavefunction(s) for diagram number 13159
    // (none)
    // Amplitude(s) for diagram number 13159
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[617], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13160 OF 15495 ***
    // Wavefunction(s) for diagram number 13160
    // (none)
    // Amplitude(s) for diagram number 13160
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[19], w_fp[214], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[18], w_fp[214], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[17], w_fp[214], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 13161 OF 15495 ***
    // Wavefunction(s) for diagram number 13161
    // (none)
    // Amplitude(s) for diagram number 13161
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[212], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13162 OF 15495 ***
    // Wavefunction(s) for diagram number 13162
    // (none)
    // Amplitude(s) for diagram number 13162
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[562], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[471], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[582], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13163 OF 15495 ***
    // Wavefunction(s) for diagram number 13163
    // (none)
    // Amplitude(s) for diagram number 13163
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[562], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[471], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[582], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];

    // *** DIAGRAM 13164 OF 15495 ***
    // Wavefunction(s) for diagram number 13164
    // (none)
    // Amplitude(s) for diagram number 13164
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[107], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[78], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[370], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 13165 OF 15495 ***
    // Wavefunction(s) for diagram number 13165
    // (none)
    // Amplitude(s) for diagram number 13165
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[4], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13166 OF 15495 ***
    // Wavefunction(s) for diagram number 13166
    // (none)
    // Amplitude(s) for diagram number 13166
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[686], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];

    // *** DIAGRAM 13167 OF 15495 ***
    // Wavefunction(s) for diagram number 13167
    // (none)
    // Amplitude(s) for diagram number 13167
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[12], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13168 OF 15495 ***
    // Wavefunction(s) for diagram number 13168
    // (none)
    // Amplitude(s) for diagram number 13168
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[686], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13169 OF 15495 ***
    // Wavefunction(s) for diagram number 13169
    // (none)
    // Amplitude(s) for diagram number 13169
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[12], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];

    // *** DIAGRAM 13170 OF 15495 ***
    // Wavefunction(s) for diagram number 13170
    // (none)
    // Amplitude(s) for diagram number 13170
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[124], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13171 OF 15495 ***
    // Wavefunction(s) for diagram number 13171
    // (none)
    // Amplitude(s) for diagram number 13171
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[261], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[352], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[278], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13172 OF 15495 ***
    // Wavefunction(s) for diagram number 13172
    // (none)
    // Amplitude(s) for diagram number 13172
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[624], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];

    // *** DIAGRAM 13173 OF 15495 ***
    // Wavefunction(s) for diagram number 13173
    // (none)
    // Amplitude(s) for diagram number 13173
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[453], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13174 OF 15495 ***
    // Wavefunction(s) for diagram number 13174
    // (none)
    // Amplitude(s) for diagram number 13174
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[214], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[214], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[214], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];

    // *** DIAGRAM 13175 OF 15495 ***
    // Wavefunction(s) for diagram number 13175
    // (none)
    // Amplitude(s) for diagram number 13175
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[4], w_fp[673], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];

    // *** DIAGRAM 13176 OF 15495 ***
    // Wavefunction(s) for diagram number 13176
    // (none)
    // Amplitude(s) for diagram number 13176
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[214], w_fp[453], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];

    // *** DIAGRAM 13177 OF 15495 ***
    // Wavefunction(s) for diagram number 13177
    // (none)
    // Amplitude(s) for diagram number 13177
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[211], w_fp[673], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13178 OF 15495 ***
    // Wavefunction(s) for diagram number 13178
    // (none)
    // Amplitude(s) for diagram number 13178
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[624], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];

    // *** DIAGRAM 13179 OF 15495 ***
    // Wavefunction(s) for diagram number 13179
    // (none)
    // Amplitude(s) for diagram number 13179
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[40], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13180 OF 15495 ***
    // Wavefunction(s) for diagram number 13180
    // (none)
    // Amplitude(s) for diagram number 13180
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[627], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13181 OF 15495 ***
    // Wavefunction(s) for diagram number 13181
    // (none)
    // Amplitude(s) for diagram number 13181
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[40], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 13182 OF 15495 ***
    // Wavefunction(s) for diagram number 13182
    // (none)
    // Amplitude(s) for diagram number 13182
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[625], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13183 OF 15495 ***
    // Wavefunction(s) for diagram number 13183
    // (none)
    // Amplitude(s) for diagram number 13183
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[214], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[214], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[214], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];

    // *** DIAGRAM 13184 OF 15495 ***
    // Wavefunction(s) for diagram number 13184
    // (none)
    // Amplitude(s) for diagram number 13184
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[214], w_fp[4], w_fp[55], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];

    // *** DIAGRAM 13185 OF 15495 ***
    // Wavefunction(s) for diagram number 13185
    // (none)
    // Amplitude(s) for diagram number 13185
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[4], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];

    // *** DIAGRAM 13186 OF 15495 ***
    // Wavefunction(s) for diagram number 13186
    // (none)
    // Amplitude(s) for diagram number 13186
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[211], w_fp[55], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13187 OF 15495 ***
    // Wavefunction(s) for diagram number 13187
    // (none)
    // Amplitude(s) for diagram number 13187
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[689], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];

    // *** DIAGRAM 13188 OF 15495 ***
    // Wavefunction(s) for diagram number 13188
    // (none)
    // Amplitude(s) for diagram number 13188
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[124], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[124], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[214], w_fp[124], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 13189 OF 15495 ***
    // Wavefunction(s) for diagram number 13189
    // (none)
    // Amplitude(s) for diagram number 13189
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[124], w_fp[745], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 13190 OF 15495 ***
    // Wavefunction(s) for diagram number 13190
    // (none)
    // Amplitude(s) for diagram number 13190
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[214], w_fp[625], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 13191 OF 15495 ***
    // Wavefunction(s) for diagram number 13191
    // (none)
    // Amplitude(s) for diagram number 13191
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[186], w_fp[689], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13192 OF 15495 ***
    // Wavefunction(s) for diagram number 13192
    // (none)
    // Amplitude(s) for diagram number 13192
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[627], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13193 OF 15495 ***
    // Wavefunction(s) for diagram number 13193
    // (none)
    // Amplitude(s) for diagram number 13193
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[214], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[126], w_fp[214], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[214], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];

    // *** DIAGRAM 13194 OF 15495 ***
    // Wavefunction(s) for diagram number 13194
    // (none)
    // Amplitude(s) for diagram number 13194
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[211], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[211], w_fp[126], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[211], w_fp[62], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13195 OF 15495 ***
    // Wavefunction(s) for diagram number 13195
    // (none)
    // Amplitude(s) for diagram number 13195
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[475], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[553], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[592], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13196 OF 15495 ***
    // Wavefunction(s) for diagram number 13196
    // (none)
    // Amplitude(s) for diagram number 13196
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[475], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[553], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[592], w_fp[1], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 13197 OF 15495 ***
    // Wavefunction(s) for diagram number 13197
    // (none)
    // Amplitude(s) for diagram number 13197
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[261], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[352], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[278], w_fp[214], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 13198 OF 15495 ***
    // Wavefunction(s) for diagram number 13198
    // (none)
    // Amplitude(s) for diagram number 13198
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[108], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[431], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[540], w_fp[280], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13199 OF 15495 ***
    // Wavefunction(s) for diagram number 13199
    // (none)
    // Amplitude(s) for diagram number 13199
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[44], w_fp[540], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[192], w_fp[540], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[75], w_fp[540], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];

    // *** DIAGRAM 13200 OF 15495 ***
    // Wavefunction(s) for diagram number 13200
    // (none)
    // Amplitude(s) for diagram number 13200
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[591], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[539], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[476], w_fp[531], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    // In CUDA, copy the local jamp to the output global-memory jamp
    constexpr int ihel0 = 0; // the allJamps buffer already points to a specific helicity _within a super-buffer for nGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcolIhelNhel( jamps, icol, ihel0, nGoodHel ) += jamp_sv[icol]; // update jamps
#else
    // In C++, copy the local jamp to the output array passed as function argument
    for( int icol = 0; icol < ncolor; icol++ )
      jamps[icol] += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
