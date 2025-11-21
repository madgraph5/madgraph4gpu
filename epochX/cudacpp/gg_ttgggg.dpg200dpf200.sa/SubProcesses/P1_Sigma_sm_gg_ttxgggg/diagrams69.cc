// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessChannelIds.h"
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
  diagramgroup69( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_RDC_DIAGRAMS
                  fptype* jamps,                  // output jamps[ncolor*2*nevt] for all events
                  const int nGoodHel,             // input: number of good helicities
#else
                  cxtype* jamps,                  // output jamps[ncolor] for this event
#endif
                  const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                  cxtype_sv* jamps,               // output jamps[ncolor*2*neppV] for this event page
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 78 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 107 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 126 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 166 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 173 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 178 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 189 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 241 );
    retrieveWf( wfs, w_cx, nevt, 242 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 247 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 250 );
    retrieveWf( wfs, w_cx, nevt, 251 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 253 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 261 );
    retrieveWf( wfs, w_cx, nevt, 264 );
    retrieveWf( wfs, w_cx, nevt, 278 );
    retrieveWf( wfs, w_cx, nevt, 280 );
    retrieveWf( wfs, w_cx, nevt, 291 );
    retrieveWf( wfs, w_cx, nevt, 305 );
    retrieveWf( wfs, w_cx, nevt, 352 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 364 );
    retrieveWf( wfs, w_cx, nevt, 365 );
    retrieveWf( wfs, w_cx, nevt, 366 );
    retrieveWf( wfs, w_cx, nevt, 370 );
    retrieveWf( wfs, w_cx, nevt, 431 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 593 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 597 );
    retrieveWf( wfs, w_cx, nevt, 598 );
    retrieveWf( wfs, w_cx, nevt, 601 );
    retrieveWf( wfs, w_cx, nevt, 602 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 605 );
    retrieveWf( wfs, w_cx, nevt, 608 );
    retrieveWf( wfs, w_cx, nevt, 611 );
    retrieveWf( wfs, w_cx, nevt, 612 );
    retrieveWf( wfs, w_cx, nevt, 615 );
    retrieveWf( wfs, w_cx, nevt, 625 );
    retrieveWf( wfs, w_cx, nevt, 630 );
    retrieveWf( wfs, w_cx, nevt, 631 );
    retrieveWf( wfs, w_cx, nevt, 632 );
    retrieveWf( wfs, w_cx, nevt, 633 );
    retrieveWf( wfs, w_cx, nevt, 634 );
    retrieveWf( wfs, w_cx, nevt, 635 );
    retrieveWf( wfs, w_cx, nevt, 639 );
    retrieveWf( wfs, w_cx, nevt, 640 );
    retrieveWf( wfs, w_cx, nevt, 641 );
    retrieveWf( wfs, w_cx, nevt, 642 );
    retrieveWf( wfs, w_cx, nevt, 643 );
    retrieveWf( wfs, w_cx, nevt, 644 );
    retrieveWf( wfs, w_cx, nevt, 658 );
    retrieveWf( wfs, w_cx, nevt, 663 );
    retrieveWf( wfs, w_cx, nevt, 673 );
    retrieveWf( wfs, w_cx, nevt, 675 );
    retrieveWf( wfs, w_cx, nevt, 696 );
    retrieveWf( wfs, w_cx, nevt, 697 );
    retrieveWf( wfs, w_cx, nevt, 698 );
    retrieveWf( wfs, w_cx, nevt, 700 );
    retrieveWf( wfs, w_cx, nevt, 724 );
    retrieveWf( wfs, w_cx, nevt, 725 );
    retrieveWf( wfs, w_cx, nevt, 734 );
    retrieveWf( wfs, w_cx, nevt, 735 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 13601 OF 15495 ***
    // Wavefunction(s) for diagram number 13601
    // (none)
    // Amplitude(s) for diagram number 13601
    FFV1_0( w_fp[256], w_fp[497], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];

    // *** DIAGRAM 13602 OF 15495 ***
    // Wavefunction(s) for diagram number 13602
    // (none)
    // Amplitude(s) for diagram number 13602
    FFV1_0( w_fp[256], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13603 OF 15495 ***
    // Wavefunction(s) for diagram number 13603
    // (none)
    // Amplitude(s) for diagram number 13603
    VVVV1_0( w_fp[516], w_fp[1], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    VVVV3_0( w_fp[516], w_fp[1], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    VVVV4_0( w_fp[516], w_fp[1], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];

    // *** DIAGRAM 13604 OF 15495 ***
    // Wavefunction(s) for diagram number 13604
    // (none)
    // Amplitude(s) for diagram number 13604
    VVV1_0( w_fp[154], w_fp[4], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];

    // *** DIAGRAM 13605 OF 15495 ***
    // Wavefunction(s) for diagram number 13605
    // (none)
    // Amplitude(s) for diagram number 13605
    VVV1_0( w_fp[1], w_fp[154], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];

    // *** DIAGRAM 13606 OF 15495 ***
    // Wavefunction(s) for diagram number 13606
    // (none)
    // Amplitude(s) for diagram number 13606
    FFV1_0( w_fp[218], w_fp[2], w_fp[700], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13607 OF 15495 ***
    // Wavefunction(s) for diagram number 13607
    // (none)
    // Amplitude(s) for diagram number 13607
    FFV1_0( w_fp[218], w_fp[497], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];

    // *** DIAGRAM 13608 OF 15495 ***
    // Wavefunction(s) for diagram number 13608
    // (none)
    // Amplitude(s) for diagram number 13608
    FFV1_0( w_fp[658], w_fp[128], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13609 OF 15495 ***
    // Wavefunction(s) for diagram number 13609
    // (none)
    // Amplitude(s) for diagram number 13609
    FFV1_0( w_fp[256], w_fp[697], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13610 OF 15495 ***
    // Wavefunction(s) for diagram number 13610
    // (none)
    // Amplitude(s) for diagram number 13610
    FFV1_0( w_fp[658], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13611 OF 15495 ***
    // Wavefunction(s) for diagram number 13611
    // (none)
    // Amplitude(s) for diagram number 13611
    FFV1_0( w_fp[256], w_fp[2], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13612 OF 15495 ***
    // Wavefunction(s) for diagram number 13612
    // (none)
    // Amplitude(s) for diagram number 13612
    VVVV1_0( w_fp[0], w_fp[361], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[361], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[361], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 13613 OF 15495 ***
    // Wavefunction(s) for diagram number 13613
    // (none)
    // Amplitude(s) for diagram number 13613
    VVV1_0( w_fp[154], w_fp[4], w_fp[238], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 13614 OF 15495 ***
    // Wavefunction(s) for diagram number 13614
    // (none)
    // Amplitude(s) for diagram number 13614
    VVV1_0( w_fp[361], w_fp[4], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 13615 OF 15495 ***
    // Wavefunction(s) for diagram number 13615
    // (none)
    // Amplitude(s) for diagram number 13615
    FFV1_0( w_fp[218], w_fp[2], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13616 OF 15495 ***
    // Wavefunction(s) for diagram number 13616
    // (none)
    // Amplitude(s) for diagram number 13616
    FFV1_0( w_fp[481], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 13617 OF 15495 ***
    // Wavefunction(s) for diagram number 13617
    // (none)
    // Amplitude(s) for diagram number 13617
    VVVV1_0( w_fp[0], w_fp[1], w_fp[154], w_fp[130], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[154], w_fp[130], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[154], w_fp[130], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13618 OF 15495 ***
    // Wavefunction(s) for diagram number 13618
    // (none)
    // Amplitude(s) for diagram number 13618
    VVV1_0( w_fp[1], w_fp[130], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13619 OF 15495 ***
    // Wavefunction(s) for diagram number 13619
    // (none)
    // Amplitude(s) for diagram number 13619
    VVV1_0( w_fp[1], w_fp[154], w_fp[525], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13620 OF 15495 ***
    // Wavefunction(s) for diagram number 13620
    // (none)
    // Amplitude(s) for diagram number 13620
    FFV1_0( w_fp[218], w_fp[697], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13621 OF 15495 ***
    // Wavefunction(s) for diagram number 13621
    // (none)
    // Amplitude(s) for diagram number 13621
    FFV1_0( w_fp[481], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13622 OF 15495 ***
    // Wavefunction(s) for diagram number 13622
    // (none)
    // Amplitude(s) for diagram number 13622
    VVV1_0( w_fp[71], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVV1_0( w_fp[132], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    VVV1_0( w_fp[65], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 13623 OF 15495 ***
    // Wavefunction(s) for diagram number 13623
    // (none)
    // Amplitude(s) for diagram number 13623
    FFV1_0( w_fp[218], w_fp[2], w_fp[71], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[132], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13624 OF 15495 ***
    // Wavefunction(s) for diagram number 13624
    // (none)
    // Amplitude(s) for diagram number 13624
    FFV1_0( w_fp[256], w_fp[2], w_fp[468], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13625 OF 15495 ***
    // Wavefunction(s) for diagram number 13625
    // (none)
    // Amplitude(s) for diagram number 13625
    VVV1_0( w_fp[468], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0( w_fp[579], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    VVV1_0( w_fp[556], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13626 OF 15495 ***
    // Wavefunction(s) for diagram number 13626
    // (none)
    // Amplitude(s) for diagram number 13626
    VVV1_0( w_fp[0], w_fp[291], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[248], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[264], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13627 OF 15495 ***
    // Wavefunction(s) for diagram number 13627
    // (none)
    // Amplitude(s) for diagram number 13627
    FFV1_0( w_fp[527], w_fp[2], w_fp[177], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[305], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13628 OF 15495 ***
    // Wavefunction(s) for diagram number 13628
    // (none)
    // Amplitude(s) for diagram number 13628
    FFV1_0( w_fp[527], w_fp[242], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    FFV1_0( w_fp[527], w_fp[47], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    FFV1_0( w_fp[527], w_fp[13], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 13629 OF 15495 ***
    // Wavefunction(s) for diagram number 13629
    // (none)
    // Amplitude(s) for diagram number 13629
    FFV1_0( w_fp[256], w_fp[2], w_fp[593], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13630 OF 15495 ***
    // Wavefunction(s) for diagram number 13630
    // (none)
    // Amplitude(s) for diagram number 13630
    VVV1_0( w_fp[593], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0( w_fp[581], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    VVV1_0( w_fp[533], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13631 OF 15495 ***
    // Wavefunction(s) for diagram number 13631
    // (none)
    // Amplitude(s) for diagram number 13631
    FFV1_0( w_fp[256], w_fp[242], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    FFV1_0( w_fp[256], w_fp[47], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    FFV1_0( w_fp[256], w_fp[13], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13632 OF 15495 ***
    // Wavefunction(s) for diagram number 13632
    // (none)
    // Amplitude(s) for diagram number 13632
    VVV1_0( w_fp[0], w_fp[177], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[305], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[10], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13633 OF 15495 ***
    // Wavefunction(s) for diagram number 13633
    // (none)
    // Amplitude(s) for diagram number 13633
    FFV1_0( w_fp[168], w_fp[2], w_fp[167], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[159], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[644], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[643], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[642], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[641], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[640], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[639], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 13634 OF 15495 ***
    // Wavefunction(s) for diagram number 13634
    // (none)
    // Amplitude(s) for diagram number 13634
    VVV1_0( w_fp[254], w_fp[7], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13635 OF 15495 ***
    // Wavefunction(s) for diagram number 13635
    // (none)
    // Amplitude(s) for diagram number 13635
    FFV1_0( w_fp[578], w_fp[2], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];

    // *** DIAGRAM 13636 OF 15495 ***
    // Wavefunction(s) for diagram number 13636
    // (none)
    // Amplitude(s) for diagram number 13636
    FFV1_0( w_fp[173], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13637 OF 15495 ***
    // Wavefunction(s) for diagram number 13637
    // (none)
    // Amplitude(s) for diagram number 13637
    FFV1_0( w_fp[578], w_fp[148], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13638 OF 15495 ***
    // Wavefunction(s) for diagram number 13638
    // (none)
    // Amplitude(s) for diagram number 13638
    FFV1_0( w_fp[173], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13639 OF 15495 ***
    // Wavefunction(s) for diagram number 13639
    // (none)
    // Amplitude(s) for diagram number 13639
    VVV1_0( w_fp[1], w_fp[69], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13640 OF 15495 ***
    // Wavefunction(s) for diagram number 13640
    // (none)
    // Amplitude(s) for diagram number 13640
    FFV1_0( w_fp[594], w_fp[2], w_fp[253], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[364], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[365], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13641 OF 15495 ***
    // Wavefunction(s) for diagram number 13641
    // (none)
    // Amplitude(s) for diagram number 13641
    FFV1_0( w_fp[255], w_fp[696], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];

    // *** DIAGRAM 13642 OF 15495 ***
    // Wavefunction(s) for diagram number 13642
    // (none)
    // Amplitude(s) for diagram number 13642
    FFV1_0( w_fp[255], w_fp[2], w_fp[598], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13643 OF 15495 ***
    // Wavefunction(s) for diagram number 13643
    // (none)
    // Amplitude(s) for diagram number 13643
    VVVV1_0( w_fp[514], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV3_0( w_fp[514], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    VVVV4_0( w_fp[514], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13644 OF 15495 ***
    // Wavefunction(s) for diagram number 13644
    // (none)
    // Amplitude(s) for diagram number 13644
    VVV1_0( w_fp[150], w_fp[7], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13645 OF 15495 ***
    // Wavefunction(s) for diagram number 13645
    // (none)
    // Amplitude(s) for diagram number 13645
    VVV1_0( w_fp[1], w_fp[150], w_fp[598], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13646 OF 15495 ***
    // Wavefunction(s) for diagram number 13646
    // (none)
    // Amplitude(s) for diagram number 13646
    FFV1_0( w_fp[176], w_fp[2], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13647 OF 15495 ***
    // Wavefunction(s) for diagram number 13647
    // (none)
    // Amplitude(s) for diagram number 13647
    FFV1_0( w_fp[176], w_fp[696], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];

    // *** DIAGRAM 13648 OF 15495 ***
    // Wavefunction(s) for diagram number 13648
    // (none)
    // Amplitude(s) for diagram number 13648
    FFV1_0( w_fp[189], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13649 OF 15495 ***
    // Wavefunction(s) for diagram number 13649
    // (none)
    // Amplitude(s) for diagram number 13649
    FFV1_0( w_fp[255], w_fp[698], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13650 OF 15495 ***
    // Wavefunction(s) for diagram number 13650
    // (none)
    // Amplitude(s) for diagram number 13650
    FFV1_0( w_fp[189], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13651 OF 15495 ***
    // Wavefunction(s) for diagram number 13651
    // (none)
    // Amplitude(s) for diagram number 13651
    FFV1_0( w_fp[255], w_fp[2], w_fp[602], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13652 OF 15495 ***
    // Wavefunction(s) for diagram number 13652
    // (none)
    // Amplitude(s) for diagram number 13652
    VVVV1_0( w_fp[0], w_fp[254], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[254], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[254], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];

    // *** DIAGRAM 13653 OF 15495 ***
    // Wavefunction(s) for diagram number 13653
    // (none)
    // Amplitude(s) for diagram number 13653
    VVV1_0( w_fp[150], w_fp[7], w_fp[182], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13654 OF 15495 ***
    // Wavefunction(s) for diagram number 13654
    // (none)
    // Amplitude(s) for diagram number 13654
    VVV1_0( w_fp[254], w_fp[7], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13655 OF 15495 ***
    // Wavefunction(s) for diagram number 13655
    // (none)
    // Amplitude(s) for diagram number 13655
    FFV1_0( w_fp[176], w_fp[2], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13656 OF 15495 ***
    // Wavefunction(s) for diagram number 13656
    // (none)
    // Amplitude(s) for diagram number 13656
    FFV1_0( w_fp[555], w_fp[2], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

    // *** DIAGRAM 13657 OF 15495 ***
    // Wavefunction(s) for diagram number 13657
    // (none)
    // Amplitude(s) for diagram number 13657
    VVVV1_0( w_fp[0], w_fp[1], w_fp[150], w_fp[69], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[150], w_fp[69], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[150], w_fp[69], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13658 OF 15495 ***
    // Wavefunction(s) for diagram number 13658
    // (none)
    // Amplitude(s) for diagram number 13658
    VVV1_0( w_fp[1], w_fp[69], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13659 OF 15495 ***
    // Wavefunction(s) for diagram number 13659
    // (none)
    // Amplitude(s) for diagram number 13659
    VVV1_0( w_fp[1], w_fp[150], w_fp[602], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13660 OF 15495 ***
    // Wavefunction(s) for diagram number 13660
    // (none)
    // Amplitude(s) for diagram number 13660
    FFV1_0( w_fp[176], w_fp[698], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13661 OF 15495 ***
    // Wavefunction(s) for diagram number 13661
    // (none)
    // Amplitude(s) for diagram number 13661
    FFV1_0( w_fp[555], w_fp[148], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13662 OF 15495 ***
    // Wavefunction(s) for diagram number 13662
    // (none)
    // Amplitude(s) for diagram number 13662
    VVV1_0( w_fp[136], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVV1_0( w_fp[137], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[551], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 13663 OF 15495 ***
    // Wavefunction(s) for diagram number 13663
    // (none)
    // Amplitude(s) for diagram number 13663
    FFV1_0( w_fp[176], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13664 OF 15495 ***
    // Wavefunction(s) for diagram number 13664
    // (none)
    // Amplitude(s) for diagram number 13664
    FFV1_0( w_fp[255], w_fp[2], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13665 OF 15495 ***
    // Wavefunction(s) for diagram number 13665
    // (none)
    // Amplitude(s) for diagram number 13665
    VVV1_0( w_fp[521], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[526], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    VVV1_0( w_fp[524], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13666 OF 15495 ***
    // Wavefunction(s) for diagram number 13666
    // (none)
    // Amplitude(s) for diagram number 13666
    VVV1_0( w_fp[0], w_fp[253], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[364], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[365], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13667 OF 15495 ***
    // Wavefunction(s) for diagram number 13667
    // (none)
    // Amplitude(s) for diagram number 13667
    VVV1_0( w_fp[363], w_fp[5], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13668 OF 15495 ***
    // Wavefunction(s) for diagram number 13668
    // (none)
    // Amplitude(s) for diagram number 13668
    FFV1_0( w_fp[577], w_fp[2], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];

    // *** DIAGRAM 13669 OF 15495 ***
    // Wavefunction(s) for diagram number 13669
    // (none)
    // Amplitude(s) for diagram number 13669
    FFV1_0( w_fp[173], w_fp[98], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13670 OF 15495 ***
    // Wavefunction(s) for diagram number 13670
    // (none)
    // Amplitude(s) for diagram number 13670
    FFV1_0( w_fp[577], w_fp[98], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13671 OF 15495 ***
    // Wavefunction(s) for diagram number 13671
    // (none)
    // Amplitude(s) for diagram number 13671
    FFV1_0( w_fp[173], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];

    // *** DIAGRAM 13672 OF 15495 ***
    // Wavefunction(s) for diagram number 13672
    // (none)
    // Amplitude(s) for diagram number 13672
    VVV1_0( w_fp[1], w_fp[103], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13673 OF 15495 ***
    // Wavefunction(s) for diagram number 13673
    // (none)
    // Amplitude(s) for diagram number 13673
    FFV1_0( w_fp[594], w_fp[2], w_fp[107], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[78], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[370], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13674 OF 15495 ***
    // Wavefunction(s) for diagram number 13674
    // (none)
    // Amplitude(s) for diagram number 13674
    FFV1_0( w_fp[255], w_fp[245], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];

    // *** DIAGRAM 13675 OF 15495 ***
    // Wavefunction(s) for diagram number 13675
    // (none)
    // Amplitude(s) for diagram number 13675
    FFV1_0( w_fp[255], w_fp[2], w_fp[612], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13676 OF 15495 ***
    // Wavefunction(s) for diagram number 13676
    // (none)
    // Amplitude(s) for diagram number 13676
    VVVV1_0( w_fp[518], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    VVVV3_0( w_fp[518], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    VVVV4_0( w_fp[518], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 13677 OF 15495 ***
    // Wavefunction(s) for diagram number 13677
    // (none)
    // Amplitude(s) for diagram number 13677
    VVV1_0( w_fp[150], w_fp[5], w_fp[611], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 13678 OF 15495 ***
    // Wavefunction(s) for diagram number 13678
    // (none)
    // Amplitude(s) for diagram number 13678
    VVV1_0( w_fp[1], w_fp[150], w_fp[612], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 13679 OF 15495 ***
    // Wavefunction(s) for diagram number 13679
    // (none)
    // Amplitude(s) for diagram number 13679
    FFV1_0( w_fp[175], w_fp[2], w_fp[611], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13680 OF 15495 ***
    // Wavefunction(s) for diagram number 13680
    // (none)
    // Amplitude(s) for diagram number 13680
    FFV1_0( w_fp[175], w_fp[245], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 13681 OF 15495 ***
    // Wavefunction(s) for diagram number 13681
    // (none)
    // Amplitude(s) for diagram number 13681
    FFV1_0( w_fp[189], w_fp[98], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13682 OF 15495 ***
    // Wavefunction(s) for diagram number 13682
    // (none)
    // Amplitude(s) for diagram number 13682
    FFV1_0( w_fp[255], w_fp[110], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13683 OF 15495 ***
    // Wavefunction(s) for diagram number 13683
    // (none)
    // Amplitude(s) for diagram number 13683
    FFV1_0( w_fp[189], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];

    // *** DIAGRAM 13684 OF 15495 ***
    // Wavefunction(s) for diagram number 13684
    // (none)
    // Amplitude(s) for diagram number 13684
    FFV1_0( w_fp[255], w_fp[2], w_fp[615], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13685 OF 15495 ***
    // Wavefunction(s) for diagram number 13685
    // (none)
    // Amplitude(s) for diagram number 13685
    VVVV1_0( w_fp[0], w_fp[363], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[363], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[363], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 13686 OF 15495 ***
    // Wavefunction(s) for diagram number 13686
    // (none)
    // Amplitude(s) for diagram number 13686
    VVV1_0( w_fp[150], w_fp[5], w_fp[61], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 13687 OF 15495 ***
    // Wavefunction(s) for diagram number 13687
    // (none)
    // Amplitude(s) for diagram number 13687
    VVV1_0( w_fp[363], w_fp[5], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 13688 OF 15495 ***
    // Wavefunction(s) for diagram number 13688
    // (none)
    // Amplitude(s) for diagram number 13688
    FFV1_0( w_fp[175], w_fp[2], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13689 OF 15495 ***
    // Wavefunction(s) for diagram number 13689
    // (none)
    // Amplitude(s) for diagram number 13689
    FFV1_0( w_fp[530], w_fp[2], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 13690 OF 15495 ***
    // Wavefunction(s) for diagram number 13690
    // (none)
    // Amplitude(s) for diagram number 13690
    VVVV1_0( w_fp[0], w_fp[1], w_fp[150], w_fp[103], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[150], w_fp[103], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[104] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[150], w_fp[103], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];

    // *** DIAGRAM 13691 OF 15495 ***
    // Wavefunction(s) for diagram number 13691
    // (none)
    // Amplitude(s) for diagram number 13691
    VVV1_0( w_fp[1], w_fp[103], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[104] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];

    // *** DIAGRAM 13692 OF 15495 ***
    // Wavefunction(s) for diagram number 13692
    // (none)
    // Amplitude(s) for diagram number 13692
    VVV1_0( w_fp[1], w_fp[150], w_fp[615], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];

    // *** DIAGRAM 13693 OF 15495 ***
    // Wavefunction(s) for diagram number 13693
    // (none)
    // Amplitude(s) for diagram number 13693
    FFV1_0( w_fp[175], w_fp[110], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13694 OF 15495 ***
    // Wavefunction(s) for diagram number 13694
    // (none)
    // Amplitude(s) for diagram number 13694
    FFV1_0( w_fp[530], w_fp[98], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13695 OF 15495 ***
    // Wavefunction(s) for diagram number 13695
    // (none)
    // Amplitude(s) for diagram number 13695
    VVV1_0( w_fp[19], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    VVV1_0( w_fp[18], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    VVV1_0( w_fp[17], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 13696 OF 15495 ***
    // Wavefunction(s) for diagram number 13696
    // (none)
    // Amplitude(s) for diagram number 13696
    FFV1_0( w_fp[175], w_fp[2], w_fp[19], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[18], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[17], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13697 OF 15495 ***
    // Wavefunction(s) for diagram number 13697
    // (none)
    // Amplitude(s) for diagram number 13697
    FFV1_0( w_fp[255], w_fp[2], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[471], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13698 OF 15495 ***
    // Wavefunction(s) for diagram number 13698
    // (none)
    // Amplitude(s) for diagram number 13698
    VVV1_0( w_fp[562], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    VVV1_0( w_fp[471], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    VVV1_0( w_fp[582], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];

    // *** DIAGRAM 13699 OF 15495 ***
    // Wavefunction(s) for diagram number 13699
    // (none)
    // Amplitude(s) for diagram number 13699
    VVV1_0( w_fp[0], w_fp[107], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[78], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[370], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 13700 OF 15495 ***
    // Wavefunction(s) for diagram number 13700
    // (none)
    // Amplitude(s) for diagram number 13700
    VVV1_0( w_fp[360], w_fp[4], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13701 OF 15495 ***
    // Wavefunction(s) for diagram number 13701
    // (none)
    // Amplitude(s) for diagram number 13701
    FFV1_0( w_fp[595], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];

    // *** DIAGRAM 13702 OF 15495 ***
    // Wavefunction(s) for diagram number 13702
    // (none)
    // Amplitude(s) for diagram number 13702
    FFV1_0( w_fp[173], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13703 OF 15495 ***
    // Wavefunction(s) for diagram number 13703
    // (none)
    // Amplitude(s) for diagram number 13703
    FFV1_0( w_fp[595], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13704 OF 15495 ***
    // Wavefunction(s) for diagram number 13704
    // (none)
    // Amplitude(s) for diagram number 13704
    FFV1_0( w_fp[173], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13705 OF 15495 ***
    // Wavefunction(s) for diagram number 13705
    // (none)
    // Amplitude(s) for diagram number 13705
    VVV1_0( w_fp[1], w_fp[124], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13706 OF 15495 ***
    // Wavefunction(s) for diagram number 13706
    // (none)
    // Amplitude(s) for diagram number 13706
    FFV1_0( w_fp[594], w_fp[2], w_fp[261], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[352], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[278], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13707 OF 15495 ***
    // Wavefunction(s) for diagram number 13707
    // (none)
    // Amplitude(s) for diagram number 13707
    FFV1_0( w_fp[255], w_fp[663], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 13708 OF 15495 ***
    // Wavefunction(s) for diagram number 13708
    // (none)
    // Amplitude(s) for diagram number 13708
    FFV1_0( w_fp[255], w_fp[2], w_fp[453], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13709 OF 15495 ***
    // Wavefunction(s) for diagram number 13709
    // (none)
    // Amplitude(s) for diagram number 13709
    VVVV1_0( w_fp[522], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    VVVV3_0( w_fp[522], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    VVVV4_0( w_fp[522], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 13710 OF 15495 ***
    // Wavefunction(s) for diagram number 13710
    // (none)
    // Amplitude(s) for diagram number 13710
    VVV1_0( w_fp[150], w_fp[4], w_fp[673], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 13711 OF 15495 ***
    // Wavefunction(s) for diagram number 13711
    // (none)
    // Amplitude(s) for diagram number 13711
    VVV1_0( w_fp[1], w_fp[150], w_fp[453], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 13712 OF 15495 ***
    // Wavefunction(s) for diagram number 13712
    // (none)
    // Amplitude(s) for diagram number 13712
    FFV1_0( w_fp[202], w_fp[2], w_fp[673], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13713 OF 15495 ***
    // Wavefunction(s) for diagram number 13713
    // (none)
    // Amplitude(s) for diagram number 13713
    FFV1_0( w_fp[202], w_fp[663], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 13714 OF 15495 ***
    // Wavefunction(s) for diagram number 13714
    // (none)
    // Amplitude(s) for diagram number 13714
    FFV1_0( w_fp[189], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13715 OF 15495 ***
    // Wavefunction(s) for diagram number 13715
    // (none)
    // Amplitude(s) for diagram number 13715
    FFV1_0( w_fp[255], w_fp[675], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13716 OF 15495 ***
    // Wavefunction(s) for diagram number 13716
    // (none)
    // Amplitude(s) for diagram number 13716
    FFV1_0( w_fp[189], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13717 OF 15495 ***
    // Wavefunction(s) for diagram number 13717
    // (none)
    // Amplitude(s) for diagram number 13717
    FFV1_0( w_fp[255], w_fp[2], w_fp[625], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13718 OF 15495 ***
    // Wavefunction(s) for diagram number 13718
    // (none)
    // Amplitude(s) for diagram number 13718
    VVVV1_0( w_fp[0], w_fp[360], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[360], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[360], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 13719 OF 15495 ***
    // Wavefunction(s) for diagram number 13719
    // (none)
    // Amplitude(s) for diagram number 13719
    VVV1_0( w_fp[150], w_fp[4], w_fp[55], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 13720 OF 15495 ***
    // Wavefunction(s) for diagram number 13720
    // (none)
    // Amplitude(s) for diagram number 13720
    VVV1_0( w_fp[360], w_fp[4], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 13721 OF 15495 ***
    // Wavefunction(s) for diagram number 13721
    // (none)
    // Amplitude(s) for diagram number 13721
    FFV1_0( w_fp[202], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13722 OF 15495 ***
    // Wavefunction(s) for diagram number 13722
    // (none)
    // Amplitude(s) for diagram number 13722
    FFV1_0( w_fp[545], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 13723 OF 15495 ***
    // Wavefunction(s) for diagram number 13723
    // (none)
    // Amplitude(s) for diagram number 13723
    VVVV1_0( w_fp[0], w_fp[1], w_fp[150], w_fp[124], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[150], w_fp[124], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[150], w_fp[124], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13724 OF 15495 ***
    // Wavefunction(s) for diagram number 13724
    // (none)
    // Amplitude(s) for diagram number 13724
    VVV1_0( w_fp[1], w_fp[124], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13725 OF 15495 ***
    // Wavefunction(s) for diagram number 13725
    // (none)
    // Amplitude(s) for diagram number 13725
    VVV1_0( w_fp[1], w_fp[150], w_fp[625], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13726 OF 15495 ***
    // Wavefunction(s) for diagram number 13726
    // (none)
    // Amplitude(s) for diagram number 13726
    FFV1_0( w_fp[202], w_fp[675], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13727 OF 15495 ***
    // Wavefunction(s) for diagram number 13727
    // (none)
    // Amplitude(s) for diagram number 13727
    FFV1_0( w_fp[545], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13728 OF 15495 ***
    // Wavefunction(s) for diagram number 13728
    // (none)
    // Amplitude(s) for diagram number 13728
    VVV1_0( w_fp[90], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVV1_0( w_fp[126], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    VVV1_0( w_fp[62], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 13729 OF 15495 ***
    // Wavefunction(s) for diagram number 13729
    // (none)
    // Amplitude(s) for diagram number 13729
    FFV1_0( w_fp[202], w_fp[2], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[202], w_fp[2], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[202], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13730 OF 15495 ***
    // Wavefunction(s) for diagram number 13730
    // (none)
    // Amplitude(s) for diagram number 13730
    FFV1_0( w_fp[255], w_fp[2], w_fp[475], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[592], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13731 OF 15495 ***
    // Wavefunction(s) for diagram number 13731
    // (none)
    // Amplitude(s) for diagram number 13731
    VVV1_0( w_fp[475], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[553], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    VVV1_0( w_fp[592], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13732 OF 15495 ***
    // Wavefunction(s) for diagram number 13732
    // (none)
    // Amplitude(s) for diagram number 13732
    VVV1_0( w_fp[0], w_fp[261], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[352], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[278], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13733 OF 15495 ***
    // Wavefunction(s) for diagram number 13733
    // (none)
    // Amplitude(s) for diagram number 13733
    FFV1_0( w_fp[594], w_fp[2], w_fp[108], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[431], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[280], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13734 OF 15495 ***
    // Wavefunction(s) for diagram number 13734
    // (none)
    // Amplitude(s) for diagram number 13734
    FFV1_0( w_fp[594], w_fp[51], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    FFV1_0( w_fp[594], w_fp[166], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    FFV1_0( w_fp[594], w_fp[241], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 13735 OF 15495 ***
    // Wavefunction(s) for diagram number 13735
    // (none)
    // Amplitude(s) for diagram number 13735
    FFV1_0( w_fp[255], w_fp[2], w_fp[591], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[539], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[531], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13736 OF 15495 ***
    // Wavefunction(s) for diagram number 13736
    // (none)
    // Amplitude(s) for diagram number 13736
    VVV1_0( w_fp[591], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[539], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    VVV1_0( w_fp[531], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13737 OF 15495 ***
    // Wavefunction(s) for diagram number 13737
    // (none)
    // Amplitude(s) for diagram number 13737
    FFV1_0( w_fp[255], w_fp[51], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    FFV1_0( w_fp[255], w_fp[166], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    FFV1_0( w_fp[255], w_fp[241], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13738 OF 15495 ***
    // Wavefunction(s) for diagram number 13738
    // (none)
    // Amplitude(s) for diagram number 13738
    VVV1_0( w_fp[0], w_fp[108], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[431], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[104] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[280], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13739 OF 15495 ***
    // Wavefunction(s) for diagram number 13739
    // (none)
    // Amplitude(s) for diagram number 13739
    FFV1_0( w_fp[174], w_fp[2], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[162], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[635], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[104] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[634], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[633], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[632], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[631], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[630], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 13740 OF 15495 ***
    // Wavefunction(s) for diagram number 13740
    // (none)
    // Amplitude(s) for diagram number 13740
    VVV1_0( w_fp[254], w_fp[6], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13741 OF 15495 ***
    // Wavefunction(s) for diagram number 13741
    // (none)
    // Amplitude(s) for diagram number 13741
    FFV1_0( w_fp[495], w_fp[2], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];

    // *** DIAGRAM 13742 OF 15495 ***
    // Wavefunction(s) for diagram number 13742
    // (none)
    // Amplitude(s) for diagram number 13742
    FFV1_0( w_fp[178], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13743 OF 15495 ***
    // Wavefunction(s) for diagram number 13743
    // (none)
    // Amplitude(s) for diagram number 13743
    FFV1_0( w_fp[495], w_fp[148], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13744 OF 15495 ***
    // Wavefunction(s) for diagram number 13744
    // (none)
    // Amplitude(s) for diagram number 13744
    FFV1_0( w_fp[178], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13745 OF 15495 ***
    // Wavefunction(s) for diagram number 13745
    // (none)
    // Amplitude(s) for diagram number 13745
    VVV1_0( w_fp[1], w_fp[68], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13746 OF 15495 ***
    // Wavefunction(s) for diagram number 13746
    // (none)
    // Amplitude(s) for diagram number 13746
    FFV1_0( w_fp[596], w_fp[2], w_fp[247], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[251], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[250], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13747 OF 15495 ***
    // Wavefunction(s) for diagram number 13747
    // (none)
    // Amplitude(s) for diagram number 13747
    FFV1_0( w_fp[252], w_fp[696], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];

    // *** DIAGRAM 13748 OF 15495 ***
    // Wavefunction(s) for diagram number 13748
    // (none)
    // Amplitude(s) for diagram number 13748
    FFV1_0( w_fp[252], w_fp[2], w_fp[597], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13749 OF 15495 ***
    // Wavefunction(s) for diagram number 13749
    // (none)
    // Amplitude(s) for diagram number 13749
    VVVV1_0( w_fp[514], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV3_0( w_fp[514], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVVV4_0( w_fp[514], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13750 OF 15495 ***
    // Wavefunction(s) for diagram number 13750
    // (none)
    // Amplitude(s) for diagram number 13750
    VVV1_0( w_fp[144], w_fp[6], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 13751 OF 15495 ***
    // Wavefunction(s) for diagram number 13751
    // (none)
    // Amplitude(s) for diagram number 13751
    VVV1_0( w_fp[1], w_fp[144], w_fp[597], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13752 OF 15495 ***
    // Wavefunction(s) for diagram number 13752
    // (none)
    // Amplitude(s) for diagram number 13752
    FFV1_0( w_fp[181], w_fp[2], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13753 OF 15495 ***
    // Wavefunction(s) for diagram number 13753
    // (none)
    // Amplitude(s) for diagram number 13753
    FFV1_0( w_fp[181], w_fp[696], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];

    // *** DIAGRAM 13754 OF 15495 ***
    // Wavefunction(s) for diagram number 13754
    // (none)
    // Amplitude(s) for diagram number 13754
    FFV1_0( w_fp[187], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13755 OF 15495 ***
    // Wavefunction(s) for diagram number 13755
    // (none)
    // Amplitude(s) for diagram number 13755
    FFV1_0( w_fp[252], w_fp[698], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13756 OF 15495 ***
    // Wavefunction(s) for diagram number 13756
    // (none)
    // Amplitude(s) for diagram number 13756
    FFV1_0( w_fp[187], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13757 OF 15495 ***
    // Wavefunction(s) for diagram number 13757
    // (none)
    // Amplitude(s) for diagram number 13757
    FFV1_0( w_fp[252], w_fp[2], w_fp[601], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13758 OF 15495 ***
    // Wavefunction(s) for diagram number 13758
    // (none)
    // Amplitude(s) for diagram number 13758
    VVVV1_0( w_fp[0], w_fp[254], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[254], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[254], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];

    // *** DIAGRAM 13759 OF 15495 ***
    // Wavefunction(s) for diagram number 13759
    // (none)
    // Amplitude(s) for diagram number 13759
    VVV1_0( w_fp[144], w_fp[6], w_fp[182], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13760 OF 15495 ***
    // Wavefunction(s) for diagram number 13760
    // (none)
    // Amplitude(s) for diagram number 13760
    VVV1_0( w_fp[254], w_fp[6], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13761 OF 15495 ***
    // Wavefunction(s) for diagram number 13761
    // (none)
    // Amplitude(s) for diagram number 13761
    FFV1_0( w_fp[181], w_fp[2], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13762 OF 15495 ***
    // Wavefunction(s) for diagram number 13762
    // (none)
    // Amplitude(s) for diagram number 13762
    FFV1_0( w_fp[583], w_fp[2], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];

    // *** DIAGRAM 13763 OF 15495 ***
    // Wavefunction(s) for diagram number 13763
    // (none)
    // Amplitude(s) for diagram number 13763
    VVVV1_0( w_fp[0], w_fp[1], w_fp[144], w_fp[68], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[144], w_fp[68], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[144], w_fp[68], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13764 OF 15495 ***
    // Wavefunction(s) for diagram number 13764
    // (none)
    // Amplitude(s) for diagram number 13764
    VVV1_0( w_fp[1], w_fp[68], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13765 OF 15495 ***
    // Wavefunction(s) for diagram number 13765
    // (none)
    // Amplitude(s) for diagram number 13765
    VVV1_0( w_fp[1], w_fp[144], w_fp[601], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13766 OF 15495 ***
    // Wavefunction(s) for diagram number 13766
    // (none)
    // Amplitude(s) for diagram number 13766
    FFV1_0( w_fp[181], w_fp[698], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13767 OF 15495 ***
    // Wavefunction(s) for diagram number 13767
    // (none)
    // Amplitude(s) for diagram number 13767
    FFV1_0( w_fp[583], w_fp[148], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13768 OF 15495 ***
    // Wavefunction(s) for diagram number 13768
    // (none)
    // Amplitude(s) for diagram number 13768
    VVV1_0( w_fp[136], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVV1_0( w_fp[137], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    VVV1_0( w_fp[551], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];

    // *** DIAGRAM 13769 OF 15495 ***
    // Wavefunction(s) for diagram number 13769
    // (none)
    // Amplitude(s) for diagram number 13769
    FFV1_0( w_fp[181], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13770 OF 15495 ***
    // Wavefunction(s) for diagram number 13770
    // (none)
    // Amplitude(s) for diagram number 13770
    FFV1_0( w_fp[252], w_fp[2], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[252], w_fp[2], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[252], w_fp[2], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13771 OF 15495 ***
    // Wavefunction(s) for diagram number 13771
    // (none)
    // Amplitude(s) for diagram number 13771
    VVV1_0( w_fp[519], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    VVV1_0( w_fp[454], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    VVV1_0( w_fp[575], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 13772 OF 15495 ***
    // Wavefunction(s) for diagram number 13772
    // (none)
    // Amplitude(s) for diagram number 13772
    VVV1_0( w_fp[0], w_fp[247], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[251], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[250], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13773 OF 15495 ***
    // Wavefunction(s) for diagram number 13773
    // (none)
    // Amplitude(s) for diagram number 13773
    VVV1_0( w_fp[362], w_fp[5], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13774 OF 15495 ***
    // Wavefunction(s) for diagram number 13774
    // (none)
    // Amplitude(s) for diagram number 13774
    FFV1_0( w_fp[438], w_fp[2], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];

    // *** DIAGRAM 13775 OF 15495 ***
    // Wavefunction(s) for diagram number 13775
    // (none)
    // Amplitude(s) for diagram number 13775
    FFV1_0( w_fp[178], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13776 OF 15495 ***
    // Wavefunction(s) for diagram number 13776
    // (none)
    // Amplitude(s) for diagram number 13776
    FFV1_0( w_fp[438], w_fp[118], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13777 OF 15495 ***
    // Wavefunction(s) for diagram number 13777
    // (none)
    // Amplitude(s) for diagram number 13777
    FFV1_0( w_fp[178], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[326] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 13778 OF 15495 ***
    // Wavefunction(s) for diagram number 13778
    // (none)
    // Amplitude(s) for diagram number 13778
    VVV1_0( w_fp[1], w_fp[87], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13779 OF 15495 ***
    // Wavefunction(s) for diagram number 13779
    // (none)
    // Amplitude(s) for diagram number 13779
    FFV1_0( w_fp[596], w_fp[2], w_fp[82], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[73], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[366], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13780 OF 15495 ***
    // Wavefunction(s) for diagram number 13780
    // (none)
    // Amplitude(s) for diagram number 13780
    FFV1_0( w_fp[252], w_fp[246], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];

    // *** DIAGRAM 13781 OF 15495 ***
    // Wavefunction(s) for diagram number 13781
    // (none)
    // Amplitude(s) for diagram number 13781
    FFV1_0( w_fp[252], w_fp[2], w_fp[605], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13782 OF 15495 ***
    // Wavefunction(s) for diagram number 13782
    // (none)
    // Amplitude(s) for diagram number 13782
    VVVV1_0( w_fp[576], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    VVVV3_0( w_fp[576], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    VVVV4_0( w_fp[576], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 13783 OF 15495 ***
    // Wavefunction(s) for diagram number 13783
    // (none)
    // Amplitude(s) for diagram number 13783
    VVV1_0( w_fp[144], w_fp[5], w_fp[604], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 13784 OF 15495 ***
    // Wavefunction(s) for diagram number 13784
    // (none)
    // Amplitude(s) for diagram number 13784
    VVV1_0( w_fp[1], w_fp[144], w_fp[605], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 13785 OF 15495 ***
    // Wavefunction(s) for diagram number 13785
    // (none)
    // Amplitude(s) for diagram number 13785
    FFV1_0( w_fp[180], w_fp[2], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13786 OF 15495 ***
    // Wavefunction(s) for diagram number 13786
    // (none)
    // Amplitude(s) for diagram number 13786
    FFV1_0( w_fp[180], w_fp[246], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 13787 OF 15495 ***
    // Wavefunction(s) for diagram number 13787
    // (none)
    // Amplitude(s) for diagram number 13787
    FFV1_0( w_fp[187], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13788 OF 15495 ***
    // Wavefunction(s) for diagram number 13788
    // (none)
    // Amplitude(s) for diagram number 13788
    FFV1_0( w_fp[252], w_fp[92], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13789 OF 15495 ***
    // Wavefunction(s) for diagram number 13789
    // (none)
    // Amplitude(s) for diagram number 13789
    FFV1_0( w_fp[187], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13790 OF 15495 ***
    // Wavefunction(s) for diagram number 13790
    // (none)
    // Amplitude(s) for diagram number 13790
    FFV1_0( w_fp[252], w_fp[2], w_fp[608], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13791 OF 15495 ***
    // Wavefunction(s) for diagram number 13791
    // (none)
    // Amplitude(s) for diagram number 13791
    VVVV1_0( w_fp[0], w_fp[362], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[362], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[362], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];

    // *** DIAGRAM 13792 OF 15495 ***
    // Wavefunction(s) for diagram number 13792
    // (none)
    // Amplitude(s) for diagram number 13792
    VVV1_0( w_fp[144], w_fp[5], w_fp[64], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 13793 OF 15495 ***
    // Wavefunction(s) for diagram number 13793
    // (none)
    // Amplitude(s) for diagram number 13793
    VVV1_0( w_fp[362], w_fp[5], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];

    // *** DIAGRAM 13794 OF 15495 ***
    // Wavefunction(s) for diagram number 13794
    // (none)
    // Amplitude(s) for diagram number 13794
    FFV1_0( w_fp[180], w_fp[2], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13795 OF 15495 ***
    // Wavefunction(s) for diagram number 13795
    // (none)
    // Amplitude(s) for diagram number 13795
    FFV1_0( w_fp[557], w_fp[2], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];

    // *** DIAGRAM 13796 OF 15495 ***
    // Wavefunction(s) for diagram number 13796
    // (none)
    // Amplitude(s) for diagram number 13796
    VVVV1_0( w_fp[0], w_fp[1], w_fp[144], w_fp[87], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[144], w_fp[87], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[144], w_fp[87], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13797 OF 15495 ***
    // Wavefunction(s) for diagram number 13797
    // (none)
    // Amplitude(s) for diagram number 13797
    VVV1_0( w_fp[1], w_fp[87], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 13798 OF 15495 ***
    // Wavefunction(s) for diagram number 13798
    // (none)
    // Amplitude(s) for diagram number 13798
    VVV1_0( w_fp[1], w_fp[144], w_fp[608], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13799 OF 15495 ***
    // Wavefunction(s) for diagram number 13799
    // (none)
    // Amplitude(s) for diagram number 13799
    FFV1_0( w_fp[180], w_fp[92], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13800 OF 15495 ***
    // Wavefunction(s) for diagram number 13800
    // (none)
    // Amplitude(s) for diagram number 13800
    FFV1_0( w_fp[557], w_fp[118], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];

#if defined MGONGPUCPP_GPUIMPL and not defined MGONGPU_RDC_DIAGRAMS
    // *** STORE JAMPS ***
    // In CUDA (DCDIAG=0), copy the local jamp to the output global-memory jamp
    constexpr int ihel0 = 0; // allJamps buffer points to a specific helicity _within a super-buffer for nGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcolIhelNhel( jamps, icol, ihel0, nGoodHel ) += jamp_sv[icol]; // update jamps
#else
    // In C++, copy the local jamp to the output array passed as function argument
    // In CUDA (DCDIAG=1), copy the local jamp to the output array passed as function argument
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
