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
  diagramgroup64( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 23 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 25 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 96 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 105 );
    retrieveWf( wfs, w_cx, nevt, 106 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 117 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 126 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 164 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 173 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 178 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 183 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 185 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 210 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 249 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 257 );
    retrieveWf( wfs, w_cx, nevt, 273 );
    retrieveWf( wfs, w_cx, nevt, 281 );
    retrieveWf( wfs, w_cx, nevt, 345 );
    retrieveWf( wfs, w_cx, nevt, 354 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 358 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 374 );
    retrieveWf( wfs, w_cx, nevt, 375 );
    retrieveWf( wfs, w_cx, nevt, 434 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 512 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 587 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 600 );
    retrieveWf( wfs, w_cx, nevt, 619 );
    retrieveWf( wfs, w_cx, nevt, 620 );
    retrieveWf( wfs, w_cx, nevt, 622 );
    retrieveWf( wfs, w_cx, nevt, 623 );
    retrieveWf( wfs, w_cx, nevt, 624 );
    retrieveWf( wfs, w_cx, nevt, 626 );
    retrieveWf( wfs, w_cx, nevt, 627 );
    retrieveWf( wfs, w_cx, nevt, 628 );
    retrieveWf( wfs, w_cx, nevt, 651 );
    retrieveWf( wfs, w_cx, nevt, 652 );
    retrieveWf( wfs, w_cx, nevt, 653 );
    retrieveWf( wfs, w_cx, nevt, 659 );
    retrieveWf( wfs, w_cx, nevt, 661 );
    retrieveWf( wfs, w_cx, nevt, 666 );
    retrieveWf( wfs, w_cx, nevt, 667 );
    retrieveWf( wfs, w_cx, nevt, 668 );
    retrieveWf( wfs, w_cx, nevt, 669 );
    retrieveWf( wfs, w_cx, nevt, 671 );
    retrieveWf( wfs, w_cx, nevt, 672 );
    retrieveWf( wfs, w_cx, nevt, 673 );
    retrieveWf( wfs, w_cx, nevt, 674 );
    retrieveWf( wfs, w_cx, nevt, 676 );
    retrieveWf( wfs, w_cx, nevt, 677 );
    retrieveWf( wfs, w_cx, nevt, 678 );
    retrieveWf( wfs, w_cx, nevt, 681 );
    retrieveWf( wfs, w_cx, nevt, 682 );
    retrieveWf( wfs, w_cx, nevt, 700 );
#endif
#endif

    // *** DIAGRAM 12601 OF 15495 ***
    // Wavefunction(s) for diagram number 12601
    // (none)
    // Amplitude(s) for diagram number 12601
    FFV1_0( w_fp[173], w_fp[163], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[350] -= amp_sv[0];

    // *** DIAGRAM 12602 OF 15495 ***
    // Wavefunction(s) for diagram number 12602
    // (none)
    // Amplitude(s) for diagram number 12602
    FFV1_0( w_fp[577], w_fp[163], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[344] -= amp_sv[0];

    // *** DIAGRAM 12603 OF 15495 ***
    // Wavefunction(s) for diagram number 12603
    // (none)
    // Amplitude(s) for diagram number 12603
    FFV1_0( w_fp[175], w_fp[659], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] -= amp_sv[0];

    // *** DIAGRAM 12604 OF 15495 ***
    // Wavefunction(s) for diagram number 12604
    // (none)
    // Amplitude(s) for diagram number 12604
    FFV1_0( w_fp[530], w_fp[512], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[282] -= amp_sv[0];

    // *** DIAGRAM 12605 OF 15495 ***
    // Wavefunction(s) for diagram number 12605
    // (none)
    // Amplitude(s) for diagram number 12605
    FFV1_0( w_fp[176], w_fp[659], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] -= amp_sv[0];

    // *** DIAGRAM 12606 OF 15495 ***
    // Wavefunction(s) for diagram number 12606
    // (none)
    // Amplitude(s) for diagram number 12606
    FFV1_0( w_fp[555], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[271] -= amp_sv[0];

    // *** DIAGRAM 12607 OF 15495 ***
    // Wavefunction(s) for diagram number 12607
    FFV1_2( w_fp[255], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[189] );
    // Amplitude(s) for diagram number 12607
    FFV1_0( w_fp[189], w_fp[158], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= amp_sv[0];

    // *** DIAGRAM 12608 OF 15495 ***
    // Wavefunction(s) for diagram number 12608
    // (none)
    // Amplitude(s) for diagram number 12608
    FFV1_0( w_fp[255], w_fp[674], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] -= amp_sv[0];

    // *** DIAGRAM 12609 OF 15495 ***
    // Wavefunction(s) for diagram number 12609
    // (none)
    // Amplitude(s) for diagram number 12609
    FFV1_0( w_fp[189], w_fp[163], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] -= amp_sv[0];

    // *** DIAGRAM 12610 OF 15495 ***
    // Wavefunction(s) for diagram number 12610
    // (none)
    // Amplitude(s) for diagram number 12610
    FFV1_0( w_fp[255], w_fp[672], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[338] -= amp_sv[0];

    // *** DIAGRAM 12611 OF 15495 ***
    // Wavefunction(s) for diagram number 12611
    // (none)
    // Amplitude(s) for diagram number 12611
    FFV1_0( w_fp[176], w_fp[674], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= amp_sv[0];

    // *** DIAGRAM 12612 OF 15495 ***
    // Wavefunction(s) for diagram number 12612
    // (none)
    // Amplitude(s) for diagram number 12612
    FFV1_0( w_fp[555], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[295] -= amp_sv[0];

    // *** DIAGRAM 12613 OF 15495 ***
    // Wavefunction(s) for diagram number 12613
    // (none)
    // Amplitude(s) for diagram number 12613
    FFV1_0( w_fp[175], w_fp[672], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];

    // *** DIAGRAM 12614 OF 15495 ***
    // Wavefunction(s) for diagram number 12614
    // (none)
    // Amplitude(s) for diagram number 12614
    FFV1_0( w_fp[530], w_fp[163], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[342] -= amp_sv[0];

    // *** DIAGRAM 12615 OF 15495 ***
    // Wavefunction(s) for diagram number 12615
    // (none)
    // Amplitude(s) for diagram number 12615
    FFV1_0( w_fp[255], w_fp[666], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12616 OF 15495 ***
    // Wavefunction(s) for diagram number 12616
    // (none)
    // Amplitude(s) for diagram number 12616
    FFV1_0( w_fp[174], w_fp[666], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];

    // *** DIAGRAM 12617 OF 15495 ***
    // Wavefunction(s) for diagram number 12617
    // (none)
    // Amplitude(s) for diagram number 12617
    FFV1_0( w_fp[109], w_fp[666], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12618 OF 15495 ***
    // Wavefunction(s) for diagram number 12618
    // (none)
    // Amplitude(s) for diagram number 12618
    FFV1_0( w_fp[594], w_fp[512], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12619 OF 15495 ***
    // Wavefunction(s) for diagram number 12619
    // (none)
    // Amplitude(s) for diagram number 12619
    FFV1_0( w_fp[594], w_fp[156], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];

    // *** DIAGRAM 12620 OF 15495 ***
    // Wavefunction(s) for diagram number 12620
    // (none)
    // Amplitude(s) for diagram number 12620
    FFV1_0( w_fp[594], w_fp[187], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12621 OF 15495 ***
    // Wavefunction(s) for diagram number 12621
    // (none)
    // Amplitude(s) for diagram number 12621
    FFV1_0( w_fp[174], w_fp[512], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];

    // *** DIAGRAM 12622 OF 15495 ***
    // Wavefunction(s) for diagram number 12622
    // (none)
    // Amplitude(s) for diagram number 12622
    FFV1_0( w_fp[255], w_fp[156], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];

    // *** DIAGRAM 12623 OF 15495 ***
    // Wavefunction(s) for diagram number 12623
    // (none)
    // Amplitude(s) for diagram number 12623
    VVV1_0( w_fp[522], w_fp[1], w_fp[178], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12624 OF 15495 ***
    // Wavefunction(s) for diagram number 12624
    // (none)
    // Amplitude(s) for diagram number 12624
    FFV1_0( w_fp[109], w_fp[512], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12625 OF 15495 ***
    // Wavefunction(s) for diagram number 12625
    // (none)
    // Amplitude(s) for diagram number 12625
    FFV1_0( w_fp[255], w_fp[187], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12626 OF 15495 ***
    // Wavefunction(s) for diagram number 12626
    // (none)
    // Amplitude(s) for diagram number 12626
    VVV1_0( w_fp[0], w_fp[360], w_fp[178], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12627 OF 15495 ***
    // Wavefunction(s) for diagram number 12627
    // (none)
    // Amplitude(s) for diagram number 12627
    FFV1_0( w_fp[174], w_fp[156], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[156], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[156], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12628 OF 15495 ***
    // Wavefunction(s) for diagram number 12628
    // (none)
    // Amplitude(s) for diagram number 12628
    FFV1_0( w_fp[252], w_fp[669], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] -= amp_sv[0];

    // *** DIAGRAM 12629 OF 15495 ***
    // Wavefunction(s) for diagram number 12629
    // (none)
    // Amplitude(s) for diagram number 12629
    FFV1_0( w_fp[252], w_fp[667], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[254] -= amp_sv[0];

    // *** DIAGRAM 12630 OF 15495 ***
    // Wavefunction(s) for diagram number 12630
    // (none)
    // Amplitude(s) for diagram number 12630
    FFV1_0( w_fp[180], w_fp[661], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] -= amp_sv[0];

    // *** DIAGRAM 12631 OF 15495 ***
    // Wavefunction(s) for diagram number 12631
    // (none)
    // Amplitude(s) for diagram number 12631
    FFV1_0( w_fp[180], w_fp[667], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= amp_sv[0];

    // *** DIAGRAM 12632 OF 15495 ***
    // Wavefunction(s) for diagram number 12632
    // (none)
    // Amplitude(s) for diagram number 12632
    FFV1_0( w_fp[181], w_fp[661], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];

    // *** DIAGRAM 12633 OF 15495 ***
    // Wavefunction(s) for diagram number 12633
    // (none)
    // Amplitude(s) for diagram number 12633
    FFV1_0( w_fp[181], w_fp[669], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] -= amp_sv[0];

    // *** DIAGRAM 12634 OF 15495 ***
    // Wavefunction(s) for diagram number 12634
    // (none)
    // Amplitude(s) for diagram number 12634
    FFV1_0( w_fp[438], w_fp[512], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[278] -= amp_sv[0];

    // *** DIAGRAM 12635 OF 15495 ***
    // Wavefunction(s) for diagram number 12635
    // (none)
    // Amplitude(s) for diagram number 12635
    FFV1_0( w_fp[495], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[272] -= amp_sv[0];

    // *** DIAGRAM 12636 OF 15495 ***
    // Wavefunction(s) for diagram number 12636
    FFV1_2( w_fp[596], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[178] );
    // Amplitude(s) for diagram number 12636
    FFV1_0( w_fp[178], w_fp[158], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] -= amp_sv[0];

    // *** DIAGRAM 12637 OF 15495 ***
    // Wavefunction(s) for diagram number 12637
    // (none)
    // Amplitude(s) for diagram number 12637
    FFV1_0( w_fp[495], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[296] -= amp_sv[0];

    // *** DIAGRAM 12638 OF 15495 ***
    // Wavefunction(s) for diagram number 12638
    // (none)
    // Amplitude(s) for diagram number 12638
    FFV1_0( w_fp[178], w_fp[161], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[326] -= amp_sv[0];

    // *** DIAGRAM 12639 OF 15495 ***
    // Wavefunction(s) for diagram number 12639
    // (none)
    // Amplitude(s) for diagram number 12639
    FFV1_0( w_fp[438], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[320] -= amp_sv[0];

    // *** DIAGRAM 12640 OF 15495 ***
    // Wavefunction(s) for diagram number 12640
    // (none)
    // Amplitude(s) for diagram number 12640
    FFV1_0( w_fp[180], w_fp[659], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] -= amp_sv[0];

    // *** DIAGRAM 12641 OF 15495 ***
    // Wavefunction(s) for diagram number 12641
    // (none)
    // Amplitude(s) for diagram number 12641
    FFV1_0( w_fp[557], w_fp[512], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[276] -= amp_sv[0];

    // *** DIAGRAM 12642 OF 15495 ***
    // Wavefunction(s) for diagram number 12642
    // (none)
    // Amplitude(s) for diagram number 12642
    FFV1_0( w_fp[181], w_fp[659], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= amp_sv[0];

    // *** DIAGRAM 12643 OF 15495 ***
    // Wavefunction(s) for diagram number 12643
    // (none)
    // Amplitude(s) for diagram number 12643
    FFV1_0( w_fp[583], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[270] -= amp_sv[0];

    // *** DIAGRAM 12644 OF 15495 ***
    // Wavefunction(s) for diagram number 12644
    FFV1_2( w_fp[252], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[187] );
    // Amplitude(s) for diagram number 12644
    FFV1_0( w_fp[187], w_fp[158], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];

    // *** DIAGRAM 12645 OF 15495 ***
    // Wavefunction(s) for diagram number 12645
    // (none)
    // Amplitude(s) for diagram number 12645
    FFV1_0( w_fp[252], w_fp[674], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] -= amp_sv[0];

    // *** DIAGRAM 12646 OF 15495 ***
    // Wavefunction(s) for diagram number 12646
    // (none)
    // Amplitude(s) for diagram number 12646
    FFV1_0( w_fp[187], w_fp[161], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= amp_sv[0];

    // *** DIAGRAM 12647 OF 15495 ***
    // Wavefunction(s) for diagram number 12647
    // (none)
    // Amplitude(s) for diagram number 12647
    FFV1_0( w_fp[252], w_fp[671], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[314] -= amp_sv[0];

    // *** DIAGRAM 12648 OF 15495 ***
    // Wavefunction(s) for diagram number 12648
    // (none)
    // Amplitude(s) for diagram number 12648
    FFV1_0( w_fp[181], w_fp[674], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= amp_sv[0];

    // *** DIAGRAM 12649 OF 15495 ***
    // Wavefunction(s) for diagram number 12649
    // (none)
    // Amplitude(s) for diagram number 12649
    FFV1_0( w_fp[583], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[294] -= amp_sv[0];

    // *** DIAGRAM 12650 OF 15495 ***
    // Wavefunction(s) for diagram number 12650
    // (none)
    // Amplitude(s) for diagram number 12650
    FFV1_0( w_fp[180], w_fp[671], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];

    // *** DIAGRAM 12651 OF 15495 ***
    // Wavefunction(s) for diagram number 12651
    // (none)
    // Amplitude(s) for diagram number 12651
    FFV1_0( w_fp[557], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[318] -= amp_sv[0];

    // *** DIAGRAM 12652 OF 15495 ***
    // Wavefunction(s) for diagram number 12652
    // (none)
    // Amplitude(s) for diagram number 12652
    FFV1_0( w_fp[252], w_fp[666], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12653 OF 15495 ***
    // Wavefunction(s) for diagram number 12653
    // (none)
    // Amplitude(s) for diagram number 12653
    FFV1_0( w_fp[179], w_fp[666], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];

    // *** DIAGRAM 12654 OF 15495 ***
    // Wavefunction(s) for diagram number 12654
    // (none)
    // Amplitude(s) for diagram number 12654
    FFV1_0( w_fp[96], w_fp[666], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12655 OF 15495 ***
    // Wavefunction(s) for diagram number 12655
    // (none)
    // Amplitude(s) for diagram number 12655
    FFV1_0( w_fp[596], w_fp[512], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12656 OF 15495 ***
    // Wavefunction(s) for diagram number 12656
    // (none)
    // Amplitude(s) for diagram number 12656
    FFV1_0( w_fp[596], w_fp[156], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[272] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];

    // *** DIAGRAM 12657 OF 15495 ***
    // Wavefunction(s) for diagram number 12657
    // (none)
    // Amplitude(s) for diagram number 12657
    FFV1_0( w_fp[596], w_fp[185], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12658 OF 15495 ***
    // Wavefunction(s) for diagram number 12658
    // (none)
    // Amplitude(s) for diagram number 12658
    FFV1_0( w_fp[179], w_fp[512], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];

    // *** DIAGRAM 12659 OF 15495 ***
    // Wavefunction(s) for diagram number 12659
    // (none)
    // Amplitude(s) for diagram number 12659
    FFV1_0( w_fp[252], w_fp[156], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];

    // *** DIAGRAM 12660 OF 15495 ***
    // Wavefunction(s) for diagram number 12660
    // (none)
    // Amplitude(s) for diagram number 12660
    VVV1_0( w_fp[559], w_fp[1], w_fp[183], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12661 OF 15495 ***
    // Wavefunction(s) for diagram number 12661
    // (none)
    // Amplitude(s) for diagram number 12661
    FFV1_0( w_fp[96], w_fp[512], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12662 OF 15495 ***
    // Wavefunction(s) for diagram number 12662
    // (none)
    // Amplitude(s) for diagram number 12662
    FFV1_0( w_fp[252], w_fp[185], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12663 OF 15495 ***
    // Wavefunction(s) for diagram number 12663
    // (none)
    // Amplitude(s) for diagram number 12663
    VVV1_0( w_fp[0], w_fp[359], w_fp[183], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12664 OF 15495 ***
    // Wavefunction(s) for diagram number 12664
    // (none)
    // Amplitude(s) for diagram number 12664
    FFV1_0( w_fp[179], w_fp[156], w_fp[106], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[156], w_fp[117], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[156], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12665 OF 15495 ***
    // Wavefunction(s) for diagram number 12665
    // (none)
    // Amplitude(s) for diagram number 12665
    VVV1_0( w_fp[359], w_fp[7], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12666 OF 15495 ***
    // Wavefunction(s) for diagram number 12666
    // (none)
    // Amplitude(s) for diagram number 12666
    FFV1_0( w_fp[3], w_fp[668], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];

    // *** DIAGRAM 12667 OF 15495 ***
    // Wavefunction(s) for diagram number 12667
    // (none)
    // Amplitude(s) for diagram number 12667
    FFV1_0( w_fp[184], w_fp[661], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12668 OF 15495 ***
    // Wavefunction(s) for diagram number 12668
    // (none)
    // Amplitude(s) for diagram number 12668
    FFV1_0( w_fp[184], w_fp[668], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12669 OF 15495 ***
    // Wavefunction(s) for diagram number 12669
    // (none)
    // Amplitude(s) for diagram number 12669
    FFV1_0( w_fp[3], w_fp[661], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 12670 OF 15495 ***
    // Wavefunction(s) for diagram number 12670
    // (none)
    // Amplitude(s) for diagram number 12670
    VVV1_0( w_fp[1], w_fp[116], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12671 OF 15495 ***
    // Wavefunction(s) for diagram number 12671
    // (none)
    // Amplitude(s) for diagram number 12671
    FFV1_0( w_fp[3], w_fp[666], w_fp[79], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[374], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[375], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12672 OF 15495 ***
    // Wavefunction(s) for diagram number 12672
    // (none)
    // Amplitude(s) for diagram number 12672
    FFV1_0( w_fp[620], w_fp[512], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];

    // *** DIAGRAM 12673 OF 15495 ***
    // Wavefunction(s) for diagram number 12673
    // (none)
    // Amplitude(s) for diagram number 12673
    FFV1_0( w_fp[3], w_fp[512], w_fp[619], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12674 OF 15495 ***
    // Wavefunction(s) for diagram number 12674
    // (none)
    // Amplitude(s) for diagram number 12674
    VVVV1_0( w_fp[559], w_fp[1], w_fp[164], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    VVVV3_0( w_fp[559], w_fp[1], w_fp[164], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    VVVV4_0( w_fp[559], w_fp[1], w_fp[164], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 12675 OF 15495 ***
    // Wavefunction(s) for diagram number 12675
    // (none)
    // Amplitude(s) for diagram number 12675
    VVV1_0( w_fp[164], w_fp[7], w_fp[120], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];

    // *** DIAGRAM 12676 OF 15495 ***
    // Wavefunction(s) for diagram number 12676
    // (none)
    // Amplitude(s) for diagram number 12676
    VVV1_0( w_fp[1], w_fp[164], w_fp[619], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 12677 OF 15495 ***
    // Wavefunction(s) for diagram number 12677
    // (none)
    // Amplitude(s) for diagram number 12677
    FFV1_0( w_fp[3], w_fp[163], w_fp[120], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12678 OF 15495 ***
    // Wavefunction(s) for diagram number 12678
    // (none)
    // Amplitude(s) for diagram number 12678
    FFV1_0( w_fp[620], w_fp[163], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];

    // *** DIAGRAM 12679 OF 15495 ***
    // Wavefunction(s) for diagram number 12679
    // (none)
    // Amplitude(s) for diagram number 12679
    FFV1_0( w_fp[184], w_fp[659], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12680 OF 15495 ***
    // Wavefunction(s) for diagram number 12680
    // (none)
    // Amplitude(s) for diagram number 12680
    FFV1_0( w_fp[623], w_fp[512], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12681 OF 15495 ***
    // Wavefunction(s) for diagram number 12681
    // (none)
    // Amplitude(s) for diagram number 12681
    FFV1_0( w_fp[3], w_fp[659], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];

    // *** DIAGRAM 12682 OF 15495 ***
    // Wavefunction(s) for diagram number 12682
    // (none)
    // Amplitude(s) for diagram number 12682
    FFV1_0( w_fp[3], w_fp[512], w_fp[622], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12683 OF 15495 ***
    // Wavefunction(s) for diagram number 12683
    // (none)
    // Amplitude(s) for diagram number 12683
    VVVV1_0( w_fp[0], w_fp[359], w_fp[164], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[359], w_fp[164], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[359], w_fp[164], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];

    // *** DIAGRAM 12684 OF 15495 ***
    // Wavefunction(s) for diagram number 12684
    // (none)
    // Amplitude(s) for diagram number 12684
    VVV1_0( w_fp[164], w_fp[7], w_fp[114], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12685 OF 15495 ***
    // Wavefunction(s) for diagram number 12685
    // (none)
    // Amplitude(s) for diagram number 12685
    VVV1_0( w_fp[359], w_fp[7], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12686 OF 15495 ***
    // Wavefunction(s) for diagram number 12686
    // (none)
    // Amplitude(s) for diagram number 12686
    FFV1_0( w_fp[3], w_fp[163], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12687 OF 15495 ***
    // Wavefunction(s) for diagram number 12687
    // (none)
    // Amplitude(s) for diagram number 12687
    FFV1_0( w_fp[3], w_fp[672], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];

    // *** DIAGRAM 12688 OF 15495 ***
    // Wavefunction(s) for diagram number 12688
    // (none)
    // Amplitude(s) for diagram number 12688
    VVVV1_0( w_fp[0], w_fp[1], w_fp[164], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[164], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[164], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 12689 OF 15495 ***
    // Wavefunction(s) for diagram number 12689
    // (none)
    // Amplitude(s) for diagram number 12689
    VVV1_0( w_fp[1], w_fp[116], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12690 OF 15495 ***
    // Wavefunction(s) for diagram number 12690
    // (none)
    // Amplitude(s) for diagram number 12690
    VVV1_0( w_fp[1], w_fp[164], w_fp[622], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 12691 OF 15495 ***
    // Wavefunction(s) for diagram number 12691
    // (none)
    // Amplitude(s) for diagram number 12691
    FFV1_0( w_fp[184], w_fp[672], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12692 OF 15495 ***
    // Wavefunction(s) for diagram number 12692
    // (none)
    // Amplitude(s) for diagram number 12692
    FFV1_0( w_fp[623], w_fp[163], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12693 OF 15495 ***
    // Wavefunction(s) for diagram number 12693
    // (none)
    // Amplitude(s) for diagram number 12693
    VVV1_0( w_fp[106], w_fp[164], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVV1_0( w_fp[117], w_fp[164], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    VVV1_0( w_fp[59], w_fp[164], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];

    // *** DIAGRAM 12694 OF 15495 ***
    // Wavefunction(s) for diagram number 12694
    // (none)
    // Amplitude(s) for diagram number 12694
    FFV1_0( w_fp[3], w_fp[163], w_fp[106], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[163], w_fp[117], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[163], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12695 OF 15495 ***
    // Wavefunction(s) for diagram number 12695
    // (none)
    // Amplitude(s) for diagram number 12695
    FFV1_0( w_fp[3], w_fp[512], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12696 OF 15495 ***
    // Wavefunction(s) for diagram number 12696
    // (none)
    // Amplitude(s) for diagram number 12696
    VVV1_0( w_fp[479], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    VVV1_0( w_fp[474], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    VVV1_0( w_fp[444], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];

    // *** DIAGRAM 12697 OF 15495 ***
    // Wavefunction(s) for diagram number 12697
    // (none)
    // Amplitude(s) for diagram number 12697
    VVV1_0( w_fp[0], w_fp[79], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[374], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[375], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12698 OF 15495 ***
    // Wavefunction(s) for diagram number 12698
    // (none)
    // Amplitude(s) for diagram number 12698
    VVV1_0( w_fp[360], w_fp[6], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12699 OF 15495 ***
    // Wavefunction(s) for diagram number 12699
    // (none)
    // Amplitude(s) for diagram number 12699
    FFV1_0( w_fp[3], w_fp[667], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];

    // *** DIAGRAM 12700 OF 15495 ***
    // Wavefunction(s) for diagram number 12700
    // (none)
    // Amplitude(s) for diagram number 12700
    FFV1_0( w_fp[186], w_fp[661], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12701 OF 15495 ***
    // Wavefunction(s) for diagram number 12701
    // (none)
    // Amplitude(s) for diagram number 12701
    FFV1_0( w_fp[186], w_fp[667], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12702 OF 15495 ***
    // Wavefunction(s) for diagram number 12702
    // (none)
    // Amplitude(s) for diagram number 12702
    FFV1_0( w_fp[3], w_fp[661], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];

    // *** DIAGRAM 12703 OF 15495 ***
    // Wavefunction(s) for diagram number 12703
    // (none)
    // Amplitude(s) for diagram number 12703
    VVV1_0( w_fp[1], w_fp[125], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12704 OF 15495 ***
    // Wavefunction(s) for diagram number 12704
    // (none)
    // Amplitude(s) for diagram number 12704
    FFV1_0( w_fp[3], w_fp[666], w_fp[358], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[81], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12705 OF 15495 ***
    // Wavefunction(s) for diagram number 12705
    // (none)
    // Amplitude(s) for diagram number 12705
    FFV1_0( w_fp[624], w_fp[512], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[276] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];

    // *** DIAGRAM 12706 OF 15495 ***
    // Wavefunction(s) for diagram number 12706
    // (none)
    // Amplitude(s) for diagram number 12706
    FFV1_0( w_fp[3], w_fp[512], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12707 OF 15495 ***
    // Wavefunction(s) for diagram number 12707
    // (none)
    // Amplitude(s) for diagram number 12707
    VVVV1_0( w_fp[522], w_fp[1], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    VVVV3_0( w_fp[522], w_fp[1], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    VVVV4_0( w_fp[522], w_fp[1], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];

    // *** DIAGRAM 12708 OF 15495 ***
    // Wavefunction(s) for diagram number 12708
    // (none)
    // Amplitude(s) for diagram number 12708
    VVV1_0( w_fp[164], w_fp[6], w_fp[673], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];

    // *** DIAGRAM 12709 OF 15495 ***
    // Wavefunction(s) for diagram number 12709
    // (none)
    // Amplitude(s) for diagram number 12709
    VVV1_0( w_fp[1], w_fp[164], w_fp[505], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];

    // *** DIAGRAM 12710 OF 15495 ***
    // Wavefunction(s) for diagram number 12710
    // (none)
    // Amplitude(s) for diagram number 12710
    FFV1_0( w_fp[3], w_fp[161], w_fp[673], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12711 OF 15495 ***
    // Wavefunction(s) for diagram number 12711
    // (none)
    // Amplitude(s) for diagram number 12711
    FFV1_0( w_fp[624], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];

    // *** DIAGRAM 12712 OF 15495 ***
    // Wavefunction(s) for diagram number 12712
    // (none)
    // Amplitude(s) for diagram number 12712
    FFV1_0( w_fp[186], w_fp[659], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12713 OF 15495 ***
    // Wavefunction(s) for diagram number 12713
    // (none)
    // Amplitude(s) for diagram number 12713
    FFV1_0( w_fp[627], w_fp[512], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12714 OF 15495 ***
    // Wavefunction(s) for diagram number 12714
    // (none)
    // Amplitude(s) for diagram number 12714
    FFV1_0( w_fp[3], w_fp[659], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];

    // *** DIAGRAM 12715 OF 15495 ***
    // Wavefunction(s) for diagram number 12715
    // (none)
    // Amplitude(s) for diagram number 12715
    FFV1_0( w_fp[3], w_fp[512], w_fp[626], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12716 OF 15495 ***
    // Wavefunction(s) for diagram number 12716
    // (none)
    // Amplitude(s) for diagram number 12716
    VVVV1_0( w_fp[0], w_fp[360], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[360], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[360], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];

    // *** DIAGRAM 12717 OF 15495 ***
    // Wavefunction(s) for diagram number 12717
    // (none)
    // Amplitude(s) for diagram number 12717
    VVV1_0( w_fp[164], w_fp[6], w_fp[55], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];

    // *** DIAGRAM 12718 OF 15495 ***
    // Wavefunction(s) for diagram number 12718
    // (none)
    // Amplitude(s) for diagram number 12718
    VVV1_0( w_fp[360], w_fp[6], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];

    // *** DIAGRAM 12719 OF 15495 ***
    // Wavefunction(s) for diagram number 12719
    // (none)
    // Amplitude(s) for diagram number 12719
    FFV1_0( w_fp[3], w_fp[161], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12720 OF 15495 ***
    // Wavefunction(s) for diagram number 12720
    // (none)
    // Amplitude(s) for diagram number 12720
    FFV1_0( w_fp[3], w_fp[671], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];

    // *** DIAGRAM 12721 OF 15495 ***
    // Wavefunction(s) for diagram number 12721
    // (none)
    // Amplitude(s) for diagram number 12721
    VVVV1_0( w_fp[0], w_fp[1], w_fp[164], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[164], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[164], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];

    // *** DIAGRAM 12722 OF 15495 ***
    // Wavefunction(s) for diagram number 12722
    // (none)
    // Amplitude(s) for diagram number 12722
    VVV1_0( w_fp[1], w_fp[125], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];

    // *** DIAGRAM 12723 OF 15495 ***
    // Wavefunction(s) for diagram number 12723
    // (none)
    // Amplitude(s) for diagram number 12723
    VVV1_0( w_fp[1], w_fp[164], w_fp[626], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];

    // *** DIAGRAM 12724 OF 15495 ***
    // Wavefunction(s) for diagram number 12724
    // (none)
    // Amplitude(s) for diagram number 12724
    FFV1_0( w_fp[186], w_fp[671], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12725 OF 15495 ***
    // Wavefunction(s) for diagram number 12725
    // (none)
    // Amplitude(s) for diagram number 12725
    FFV1_0( w_fp[627], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12726 OF 15495 ***
    // Wavefunction(s) for diagram number 12726
    // (none)
    // Amplitude(s) for diagram number 12726
    VVV1_0( w_fp[90], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    VVV1_0( w_fp[126], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    VVV1_0( w_fp[62], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];

    // *** DIAGRAM 12727 OF 15495 ***
    // Wavefunction(s) for diagram number 12727
    // (none)
    // Amplitude(s) for diagram number 12727
    FFV1_0( w_fp[3], w_fp[161], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[161], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[161], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12728 OF 15495 ***
    // Wavefunction(s) for diagram number 12728
    // (none)
    // Amplitude(s) for diagram number 12728
    FFV1_0( w_fp[3], w_fp[512], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12729 OF 15495 ***
    // Wavefunction(s) for diagram number 12729
    // (none)
    // Amplitude(s) for diagram number 12729
    VVV1_0( w_fp[523], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    VVV1_0( w_fp[580], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    VVV1_0( w_fp[528], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];

    // *** DIAGRAM 12730 OF 15495 ***
    // Wavefunction(s) for diagram number 12730
    // (none)
    // Amplitude(s) for diagram number 12730
    VVV1_0( w_fp[0], w_fp[358], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[81], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[46], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];

    // *** DIAGRAM 12731 OF 15495 ***
    // Wavefunction(s) for diagram number 12731
    // (none)
    // Amplitude(s) for diagram number 12731
    VVV1_0( w_fp[361], w_fp[5], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12732 OF 15495 ***
    // Wavefunction(s) for diagram number 12732
    // (none)
    // Amplitude(s) for diagram number 12732
    FFV1_0( w_fp[3], w_fp[669], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];

    // *** DIAGRAM 12733 OF 15495 ***
    // Wavefunction(s) for diagram number 12733
    // (none)
    // Amplitude(s) for diagram number 12733
    FFV1_0( w_fp[188], w_fp[661], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12734 OF 15495 ***
    // Wavefunction(s) for diagram number 12734
    // (none)
    // Amplitude(s) for diagram number 12734
    FFV1_0( w_fp[188], w_fp[669], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12735 OF 15495 ***
    // Wavefunction(s) for diagram number 12735
    // (none)
    // Amplitude(s) for diagram number 12735
    FFV1_0( w_fp[3], w_fp[661], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 12736 OF 15495 ***
    // Wavefunction(s) for diagram number 12736
    // (none)
    // Amplitude(s) for diagram number 12736
    VVV1_0( w_fp[1], w_fp[131], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12737 OF 15495 ***
    // Wavefunction(s) for diagram number 12737
    // (none)
    // Amplitude(s) for diagram number 12737
    FFV1_0( w_fp[3], w_fp[666], w_fp[257], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[249], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[354], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12738 OF 15495 ***
    // Wavefunction(s) for diagram number 12738
    // (none)
    // Amplitude(s) for diagram number 12738
    FFV1_0( w_fp[442], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[270] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];

    // *** DIAGRAM 12739 OF 15495 ***
    // Wavefunction(s) for diagram number 12739
    // (none)
    // Amplitude(s) for diagram number 12739
    FFV1_0( w_fp[3], w_fp[512], w_fp[587], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12740 OF 15495 ***
    // Wavefunction(s) for diagram number 12740
    // (none)
    // Amplitude(s) for diagram number 12740
    VVVV1_0( w_fp[516], w_fp[1], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    VVVV3_0( w_fp[516], w_fp[1], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[256] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    VVVV4_0( w_fp[516], w_fp[1], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];

    // *** DIAGRAM 12741 OF 15495 ***
    // Wavefunction(s) for diagram number 12741
    // (none)
    // Amplitude(s) for diagram number 12741
    VVV1_0( w_fp[164], w_fp[5], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];

    // *** DIAGRAM 12742 OF 15495 ***
    // Wavefunction(s) for diagram number 12742
    // (none)
    // Amplitude(s) for diagram number 12742
    VVV1_0( w_fp[1], w_fp[164], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];

    // *** DIAGRAM 12743 OF 15495 ***
    // Wavefunction(s) for diagram number 12743
    // (none)
    // Amplitude(s) for diagram number 12743
    FFV1_0( w_fp[3], w_fp[158], w_fp[700], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12744 OF 15495 ***
    // Wavefunction(s) for diagram number 12744
    // (none)
    // Amplitude(s) for diagram number 12744
    FFV1_0( w_fp[442], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];

    // *** DIAGRAM 12745 OF 15495 ***
    // Wavefunction(s) for diagram number 12745
    // (none)
    // Amplitude(s) for diagram number 12745
    FFV1_0( w_fp[188], w_fp[659], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12746 OF 15495 ***
    // Wavefunction(s) for diagram number 12746
    // (none)
    // Amplitude(s) for diagram number 12746
    FFV1_0( w_fp[628], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12747 OF 15495 ***
    // Wavefunction(s) for diagram number 12747
    // (none)
    // Amplitude(s) for diagram number 12747
    FFV1_0( w_fp[3], w_fp[659], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];

    // *** DIAGRAM 12748 OF 15495 ***
    // Wavefunction(s) for diagram number 12748
    // (none)
    // Amplitude(s) for diagram number 12748
    FFV1_0( w_fp[3], w_fp[512], w_fp[600], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12749 OF 15495 ***
    // Wavefunction(s) for diagram number 12749
    // (none)
    // Amplitude(s) for diagram number 12749
    VVVV1_0( w_fp[0], w_fp[361], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[361], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[361], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];

    // *** DIAGRAM 12750 OF 15495 ***
    // Wavefunction(s) for diagram number 12750
    // (none)
    // Amplitude(s) for diagram number 12750
    VVV1_0( w_fp[164], w_fp[5], w_fp[238], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];

    // *** DIAGRAM 12751 OF 15495 ***
    // Wavefunction(s) for diagram number 12751
    // (none)
    // Amplitude(s) for diagram number 12751
    VVV1_0( w_fp[361], w_fp[5], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];

    // *** DIAGRAM 12752 OF 15495 ***
    // Wavefunction(s) for diagram number 12752
    // (none)
    // Amplitude(s) for diagram number 12752
    FFV1_0( w_fp[3], w_fp[158], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12753 OF 15495 ***
    // Wavefunction(s) for diagram number 12753
    // (none)
    // Amplitude(s) for diagram number 12753
    FFV1_0( w_fp[3], w_fp[674], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];

    // *** DIAGRAM 12754 OF 15495 ***
    // Wavefunction(s) for diagram number 12754
    // (none)
    // Amplitude(s) for diagram number 12754
    VVVV1_0( w_fp[0], w_fp[1], w_fp[164], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[164], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[164], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 12755 OF 15495 ***
    // Wavefunction(s) for diagram number 12755
    // (none)
    // Amplitude(s) for diagram number 12755
    VVV1_0( w_fp[1], w_fp[131], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12756 OF 15495 ***
    // Wavefunction(s) for diagram number 12756
    // (none)
    // Amplitude(s) for diagram number 12756
    VVV1_0( w_fp[1], w_fp[164], w_fp[600], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 12757 OF 15495 ***
    // Wavefunction(s) for diagram number 12757
    // (none)
    // Amplitude(s) for diagram number 12757
    FFV1_0( w_fp[188], w_fp[674], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12758 OF 15495 ***
    // Wavefunction(s) for diagram number 12758
    // (none)
    // Amplitude(s) for diagram number 12758
    FFV1_0( w_fp[628], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12759 OF 15495 ***
    // Wavefunction(s) for diagram number 12759
    // (none)
    // Amplitude(s) for diagram number 12759
    VVV1_0( w_fp[71], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    VVV1_0( w_fp[132], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    VVV1_0( w_fp[65], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];

    // *** DIAGRAM 12760 OF 15495 ***
    // Wavefunction(s) for diagram number 12760
    // (none)
    // Amplitude(s) for diagram number 12760
    FFV1_0( w_fp[3], w_fp[158], w_fp[71], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[158], w_fp[132], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[158], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12761 OF 15495 ***
    // Wavefunction(s) for diagram number 12761
    // (none)
    // Amplitude(s) for diagram number 12761
    FFV1_0( w_fp[3], w_fp[512], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12762 OF 15495 ***
    // Wavefunction(s) for diagram number 12762
    // (none)
    // Amplitude(s) for diagram number 12762
    VVV1_0( w_fp[549], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    VVV1_0( w_fp[451], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    VVV1_0( w_fp[571], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];

    // *** DIAGRAM 12763 OF 15495 ***
    // Wavefunction(s) for diagram number 12763
    // (none)
    // Amplitude(s) for diagram number 12763
    VVV1_0( w_fp[0], w_fp[257], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[249], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[354], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12764 OF 15495 ***
    // Wavefunction(s) for diagram number 12764
    // (none)
    // Amplitude(s) for diagram number 12764
    FFV1_0( w_fp[3], w_fp[666], w_fp[345], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[95], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[434], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12765 OF 15495 ***
    // Wavefunction(s) for diagram number 12765
    // (none)
    // Amplitude(s) for diagram number 12765
    FFV1_0( w_fp[26], w_fp[666], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    FFV1_0( w_fp[160], w_fp[666], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    FFV1_0( w_fp[105], w_fp[666], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];

    // *** DIAGRAM 12766 OF 15495 ***
    // Wavefunction(s) for diagram number 12766
    // (none)
    // Amplitude(s) for diagram number 12766
    FFV1_0( w_fp[3], w_fp[512], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12767 OF 15495 ***
    // Wavefunction(s) for diagram number 12767
    // (none)
    // Amplitude(s) for diagram number 12767
    VVV1_0( w_fp[452], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    VVV1_0( w_fp[488], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    VVV1_0( w_fp[437], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];

    // *** DIAGRAM 12768 OF 15495 ***
    // Wavefunction(s) for diagram number 12768
    // (none)
    // Amplitude(s) for diagram number 12768
    FFV1_0( w_fp[26], w_fp[512], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    FFV1_0( w_fp[160], w_fp[512], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    FFV1_0( w_fp[105], w_fp[512], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];

    // *** DIAGRAM 12769 OF 15495 ***
    // Wavefunction(s) for diagram number 12769
    // (none)
    // Amplitude(s) for diagram number 12769
    VVV1_0( w_fp[0], w_fp[345], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[95], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[434], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];

    // *** DIAGRAM 12770 OF 15495 ***
    // Wavefunction(s) for diagram number 12770
    // (none)
    // Amplitude(s) for diagram number 12770
    FFV1_0( w_fp[3], w_fp[156], w_fp[23], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[24], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[25], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[653], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[652], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[651], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];

    // *** DIAGRAM 12771 OF 15495 ***
    // Wavefunction(s) for diagram number 12771
    // (none)
    // Amplitude(s) for diagram number 12771
    FFV1_0( w_fp[355], w_fp[677], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= amp_sv[0];

    // *** DIAGRAM 12772 OF 15495 ***
    // Wavefunction(s) for diagram number 12772
    // (none)
    // Amplitude(s) for diagram number 12772
    FFV1_0( w_fp[355], w_fp[678], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] -= amp_sv[0];

    // *** DIAGRAM 12773 OF 15495 ***
    // Wavefunction(s) for diagram number 12773
    FFV1_1( w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[105] );
    // Amplitude(s) for diagram number 12773
    FFV1_0( w_fp[198], w_fp[105], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] -= amp_sv[0];

    // *** DIAGRAM 12774 OF 15495 ***
    // Wavefunction(s) for diagram number 12774
    // (none)
    // Amplitude(s) for diagram number 12774
    FFV1_0( w_fp[198], w_fp[678], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] -= amp_sv[0];

    // *** DIAGRAM 12775 OF 15495 ***
    // Wavefunction(s) for diagram number 12775
    // (none)
    // Amplitude(s) for diagram number 12775
    FFV1_0( w_fp[199], w_fp[105], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= amp_sv[0];

    // *** DIAGRAM 12776 OF 15495 ***
    // Wavefunction(s) for diagram number 12776
    // (none)
    // Amplitude(s) for diagram number 12776
    FFV1_0( w_fp[199], w_fp[677], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] -= amp_sv[0];

    // *** DIAGRAM 12777 OF 15495 ***
    // Wavefunction(s) for diagram number 12777
    // (none)
    // Amplitude(s) for diagram number 12777
    FFV1_0( w_fp[560], w_fp[477], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[406] -= amp_sv[0];

    // *** DIAGRAM 12778 OF 15495 ***
    // Wavefunction(s) for diagram number 12778
    // (none)
    // Amplitude(s) for diagram number 12778
    FFV1_0( w_fp[45], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[400] -= amp_sv[0];

    // *** DIAGRAM 12779 OF 15495 ***
    // Wavefunction(s) for diagram number 12779
    FFV1_2( w_fp[536], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[160] );
    // Amplitude(s) for diagram number 12779
    FFV1_0( w_fp[160], w_fp[191], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] -= amp_sv[0];

    // *** DIAGRAM 12780 OF 15495 ***
    // Wavefunction(s) for diagram number 12780
    // (none)
    // Amplitude(s) for diagram number 12780
    FFV1_0( w_fp[45], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[442] -= amp_sv[0];

    // *** DIAGRAM 12781 OF 15495 ***
    // Wavefunction(s) for diagram number 12781
    // (none)
    // Amplitude(s) for diagram number 12781
    FFV1_0( w_fp[160], w_fp[193], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] -= amp_sv[0];

    // *** DIAGRAM 12782 OF 15495 ***
    // Wavefunction(s) for diagram number 12782
    // (none)
    // Amplitude(s) for diagram number 12782
    FFV1_0( w_fp[560], w_fp[193], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[466] -= amp_sv[0];

    // *** DIAGRAM 12783 OF 15495 ***
    // Wavefunction(s) for diagram number 12783
    FFV1_1( w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[26] );
    // Amplitude(s) for diagram number 12783
    FFV1_0( w_fp[198], w_fp[26], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[389] -= amp_sv[0];

    // *** DIAGRAM 12784 OF 15495 ***
    // Wavefunction(s) for diagram number 12784
    // (none)
    // Amplitude(s) for diagram number 12784
    FFV1_0( w_fp[541], w_fp[477], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[403] -= amp_sv[0];

    // *** DIAGRAM 12785 OF 15495 ***
    // Wavefunction(s) for diagram number 12785
    // (none)
    // Amplitude(s) for diagram number 12785
    FFV1_0( w_fp[199], w_fp[26], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] -= amp_sv[0];

    // *** DIAGRAM 12786 OF 15495 ***
    // Wavefunction(s) for diagram number 12786
    // (none)
    // Amplitude(s) for diagram number 12786
    FFV1_0( w_fp[532], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[397] -= amp_sv[0];

    // *** DIAGRAM 12787 OF 15495 ***
    // Wavefunction(s) for diagram number 12787
    FFV1_2( w_fp[355], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[674] );
    // Amplitude(s) for diagram number 12787
    FFV1_0( w_fp[674], w_fp[191], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= amp_sv[0];

    // *** DIAGRAM 12788 OF 15495 ***
    // Wavefunction(s) for diagram number 12788
    // (none)
    // Amplitude(s) for diagram number 12788
    FFV1_0( w_fp[355], w_fp[681], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] -= amp_sv[0];

    // *** DIAGRAM 12789 OF 15495 ***
    // Wavefunction(s) for diagram number 12789
    // (none)
    // Amplitude(s) for diagram number 12789
    FFV1_0( w_fp[674], w_fp[193], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] -= amp_sv[0];

    // *** DIAGRAM 12790 OF 15495 ***
    // Wavefunction(s) for diagram number 12790
    // (none)
    // Amplitude(s) for diagram number 12790
    FFV1_0( w_fp[355], w_fp[682], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] -= amp_sv[0];

    // *** DIAGRAM 12791 OF 15495 ***
    // Wavefunction(s) for diagram number 12791
    // (none)
    // Amplitude(s) for diagram number 12791
    FFV1_0( w_fp[199], w_fp[681], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= amp_sv[0];

    // *** DIAGRAM 12792 OF 15495 ***
    // Wavefunction(s) for diagram number 12792
    // (none)
    // Amplitude(s) for diagram number 12792
    FFV1_0( w_fp[532], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[439] -= amp_sv[0];

    // *** DIAGRAM 12793 OF 15495 ***
    // Wavefunction(s) for diagram number 12793
    // (none)
    // Amplitude(s) for diagram number 12793
    FFV1_0( w_fp[198], w_fp[682], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] -= amp_sv[0];

    // *** DIAGRAM 12794 OF 15495 ***
    // Wavefunction(s) for diagram number 12794
    // (none)
    // Amplitude(s) for diagram number 12794
    FFV1_0( w_fp[541], w_fp[193], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[463] -= amp_sv[0];

    // *** DIAGRAM 12795 OF 15495 ***
    // Wavefunction(s) for diagram number 12795
    // (none)
    // Amplitude(s) for diagram number 12795
    FFV1_0( w_fp[355], w_fp[676], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12796 OF 15495 ***
    // Wavefunction(s) for diagram number 12796
    // (none)
    // Amplitude(s) for diagram number 12796
    FFV1_0( w_fp[196], w_fp[676], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];

    // *** DIAGRAM 12797 OF 15495 ***
    // Wavefunction(s) for diagram number 12797
    // (none)
    // Amplitude(s) for diagram number 12797
    FFV1_0( w_fp[129], w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12798 OF 15495 ***
    // Wavefunction(s) for diagram number 12798
    // (none)
    // Amplitude(s) for diagram number 12798
    FFV1_0( w_fp[536], w_fp[477], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12799 OF 15495 ***
    // Wavefunction(s) for diagram number 12799
    // (none)
    // Amplitude(s) for diagram number 12799
    FFV1_0( w_fp[536], w_fp[169], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];

    // *** DIAGRAM 12800 OF 15495 ***
    // Wavefunction(s) for diagram number 12800
    // (none)
    // Amplitude(s) for diagram number 12800
    FFV1_0( w_fp[536], w_fp[210], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 26 );
    storeWf( wfs, w_cx, nevt, 105 );
    storeWf( wfs, w_cx, nevt, 160 );
    storeWf( wfs, w_cx, nevt, 178 );
    storeWf( wfs, w_cx, nevt, 187 );
    storeWf( wfs, w_cx, nevt, 189 );
    storeWf( wfs, w_cx, nevt, 674 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
