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
  diagramgroup24( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 112 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 224 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 231 );
    retrieveWf( wfs, w_cx, nevt, 232 );
    retrieveWf( wfs, w_cx, nevt, 235 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 448 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 587 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 589 );
    retrieveWf( wfs, w_cx, nevt, 590 );
#endif
#endif

    // *** DIAGRAM 4601 OF 15495 ***
    // Wavefunction(s) for diagram number 4601
    // (none)
    // Amplitude(s) for diagram number 4601
    FFV1_0( w_fp[3], w_fp[224], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4602 OF 15495 ***
    // Wavefunction(s) for diagram number 4602
    // (none)
    // Amplitude(s) for diagram number 4602
    FFV1_0( w_fp[186], w_fp[197], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4603 OF 15495 ***
    // Wavefunction(s) for diagram number 4603
    // (none)
    // Amplitude(s) for diagram number 4603
    FFV1_0( w_fp[3], w_fp[541], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4604 OF 15495 ***
    // Wavefunction(s) for diagram number 4604
    // (none)
    // Amplitude(s) for diagram number 4604
    FFV1_0( w_fp[186], w_fp[541], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 4605 OF 15495 ***
    // Wavefunction(s) for diagram number 4605
    // (none)
    // Amplitude(s) for diagram number 4605
    FFV1_0( w_fp[479], w_fp[476], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[520] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];

    // *** DIAGRAM 4606 OF 15495 ***
    // Wavefunction(s) for diagram number 4606
    // (none)
    // Amplitude(s) for diagram number 4606
    FFV1_0( w_fp[479], w_fp[197], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4607 OF 15495 ***
    // Wavefunction(s) for diagram number 4607
    // (none)
    // Amplitude(s) for diagram number 4607
    FFV1_0( w_fp[479], w_fp[224], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];

    // *** DIAGRAM 4608 OF 15495 ***
    // Wavefunction(s) for diagram number 4608
    // (none)
    // Amplitude(s) for diagram number 4608
    FFV1_0( w_fp[3], w_fp[476], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4609 OF 15495 ***
    // Wavefunction(s) for diagram number 4609
    // (none)
    // Amplitude(s) for diagram number 4609
    VVV1_0( w_fp[474], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 4610 OF 15495 ***
    // Wavefunction(s) for diagram number 4610
    // (none)
    // Amplitude(s) for diagram number 4610
    FFV1_0( w_fp[186], w_fp[476], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];

    // *** DIAGRAM 4611 OF 15495 ***
    // Wavefunction(s) for diagram number 4611
    // (none)
    // Amplitude(s) for diagram number 4611
    VVV1_0( w_fp[435], w_fp[360], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 4612 OF 15495 ***
    // Wavefunction(s) for diagram number 4612
    // (none)
    // Amplitude(s) for diagram number 4612
    FFV1_0( w_fp[3], w_fp[197], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 4613 OF 15495 ***
    // Wavefunction(s) for diagram number 4613
    // (none)
    // Amplitude(s) for diagram number 4613
    VVVV1_0( w_fp[546], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    VVVV3_0( w_fp[546], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    VVVV4_0( w_fp[546], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4614 OF 15495 ***
    // Wavefunction(s) for diagram number 4614
    // (none)
    // Amplitude(s) for diagram number 4614
    VVV1_0( w_fp[228], w_fp[6], w_fp[549], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 4615 OF 15495 ***
    // Wavefunction(s) for diagram number 4615
    // (none)
    // Amplitude(s) for diagram number 4615
    VVV1_0( w_fp[228], w_fp[5], w_fp[561], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4616 OF 15495 ***
    // Wavefunction(s) for diagram number 4616
    // (none)
    // Amplitude(s) for diagram number 4616
    FFV1_0( w_fp[529], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[691] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];

    // *** DIAGRAM 4617 OF 15495 ***
    // Wavefunction(s) for diagram number 4617
    // (none)
    // Amplitude(s) for diagram number 4617
    FFV1_0( w_fp[3], w_fp[226], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4618 OF 15495 ***
    // Wavefunction(s) for diagram number 4618
    // (none)
    // Amplitude(s) for diagram number 4618
    FFV1_0( w_fp[529], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[715] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 4619 OF 15495 ***
    // Wavefunction(s) for diagram number 4619
    // (none)
    // Amplitude(s) for diagram number 4619
    FFV1_0( w_fp[3], w_fp[227], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4620 OF 15495 ***
    // Wavefunction(s) for diagram number 4620
    // (none)
    // Amplitude(s) for diagram number 4620
    FFV1_0( w_fp[447], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4621 OF 15495 ***
    // Wavefunction(s) for diagram number 4621
    // (none)
    // Amplitude(s) for diagram number 4621
    FFV1_0( w_fp[505], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4622 OF 15495 ***
    // Wavefunction(s) for diagram number 4622
    // (none)
    // Amplitude(s) for diagram number 4622
    FFV1_0( w_fp[528], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4623 OF 15495 ***
    // Wavefunction(s) for diagram number 4623
    // (none)
    // Amplitude(s) for diagram number 4623
    FFV1_0( w_fp[505], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4624 OF 15495 ***
    // Wavefunction(s) for diagram number 4624
    // (none)
    // Amplitude(s) for diagram number 4624
    FFV1_0( w_fp[528], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4625 OF 15495 ***
    // Wavefunction(s) for diagram number 4625
    // (none)
    // Amplitude(s) for diagram number 4625
    FFV1_0( w_fp[447], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4626 OF 15495 ***
    // Wavefunction(s) for diagram number 4626
    // (none)
    // Amplitude(s) for diagram number 4626
    FFV1_0( w_fp[495], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];

    // *** DIAGRAM 4627 OF 15495 ***
    // Wavefunction(s) for diagram number 4627
    // (none)
    // Amplitude(s) for diagram number 4627
    FFV1_0( w_fp[3], w_fp[473], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4628 OF 15495 ***
    // Wavefunction(s) for diagram number 4628
    // (none)
    // Amplitude(s) for diagram number 4628
    VVVV1_0( w_fp[451], w_fp[1], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    VVVV3_0( w_fp[451], w_fp[1], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV4_0( w_fp[451], w_fp[1], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 4629 OF 15495 ***
    // Wavefunction(s) for diagram number 4629
    // (none)
    // Amplitude(s) for diagram number 4629
    VVV1_0( w_fp[228], w_fp[6], w_fp[553], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 4630 OF 15495 ***
    // Wavefunction(s) for diagram number 4630
    // (none)
    // Amplitude(s) for diagram number 4630
    VVV1_0( w_fp[1], w_fp[228], w_fp[514], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 4631 OF 15495 ***
    // Wavefunction(s) for diagram number 4631
    // (none)
    // Amplitude(s) for diagram number 4631
    FFV1_0( w_fp[3], w_fp[227], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4632 OF 15495 ***
    // Wavefunction(s) for diagram number 4632
    // (none)
    // Amplitude(s) for diagram number 4632
    FFV1_0( w_fp[495], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 4633 OF 15495 ***
    // Wavefunction(s) for diagram number 4633
    // (none)
    // Amplitude(s) for diagram number 4633
    FFV1_0( w_fp[450], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];

    // *** DIAGRAM 4634 OF 15495 ***
    // Wavefunction(s) for diagram number 4634
    // (none)
    // Amplitude(s) for diagram number 4634
    FFV1_0( w_fp[3], w_fp[473], w_fp[497], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4635 OF 15495 ***
    // Wavefunction(s) for diagram number 4635
    // (none)
    // Amplitude(s) for diagram number 4635
    VVVV1_0( w_fp[488], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 4636 OF 15495 ***
    // Wavefunction(s) for diagram number 4636
    // (none)
    // Amplitude(s) for diagram number 4636
    VVV1_0( w_fp[228], w_fp[5], w_fp[571], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4637 OF 15495 ***
    // Wavefunction(s) for diagram number 4637
    // (none)
    // Amplitude(s) for diagram number 4637
    VVV1_0( w_fp[1], w_fp[228], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 4638 OF 15495 ***
    // Wavefunction(s) for diagram number 4638
    // (none)
    // Amplitude(s) for diagram number 4638
    FFV1_0( w_fp[3], w_fp[226], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4639 OF 15495 ***
    // Wavefunction(s) for diagram number 4639
    // (none)
    // Amplitude(s) for diagram number 4639
    FFV1_0( w_fp[450], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 4640 OF 15495 ***
    // Wavefunction(s) for diagram number 4640
    // (none)
    // Amplitude(s) for diagram number 4640
    VVV1_0( w_fp[576], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0( w_fp[577], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[578], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 4641 OF 15495 ***
    // Wavefunction(s) for diagram number 4641
    // (none)
    // Amplitude(s) for diagram number 4641
    FFV1_0( w_fp[3], w_fp[227], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4642 OF 15495 ***
    // Wavefunction(s) for diagram number 4642
    // (none)
    // Amplitude(s) for diagram number 4642
    VVV1_0( w_fp[579], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVV1_0( w_fp[580], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVV1_0( w_fp[581], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4643 OF 15495 ***
    // Wavefunction(s) for diagram number 4643
    // (none)
    // Amplitude(s) for diagram number 4643
    FFV1_0( w_fp[3], w_fp[226], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4644 OF 15495 ***
    // Wavefunction(s) for diagram number 4644
    // (none)
    // Amplitude(s) for diagram number 4644
    FFV1_0( w_fp[3], w_fp[473], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4645 OF 15495 ***
    // Wavefunction(s) for diagram number 4645
    // (none)
    // Amplitude(s) for diagram number 4645
    VVV1_0( w_fp[515], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[454], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    VVV1_0( w_fp[438], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 4646 OF 15495 ***
    // Wavefunction(s) for diagram number 4646
    // (none)
    // Amplitude(s) for diagram number 4646
    VVV1_0( w_fp[546], w_fp[231], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4647 OF 15495 ***
    // Wavefunction(s) for diagram number 4647
    // (none)
    // Amplitude(s) for diagram number 4647
    FFV1_0( w_fp[168], w_fp[227], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4648 OF 15495 ***
    // Wavefunction(s) for diagram number 4648
    // (none)
    // Amplitude(s) for diagram number 4648
    FFV1_0( w_fp[170], w_fp[215], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];

    // *** DIAGRAM 4649 OF 15495 ***
    // Wavefunction(s) for diagram number 4649
    // (none)
    // Amplitude(s) for diagram number 4649
    FFV1_0( w_fp[256], w_fp[543], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4650 OF 15495 ***
    // Wavefunction(s) for diagram number 4650
    // (none)
    // Amplitude(s) for diagram number 4650
    FFV1_0( w_fp[170], w_fp[543], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4651 OF 15495 ***
    // Wavefunction(s) for diagram number 4651
    // (none)
    // Amplitude(s) for diagram number 4651
    FFV1_0( w_fp[449], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4652 OF 15495 ***
    // Wavefunction(s) for diagram number 4652
    // (none)
    // Amplitude(s) for diagram number 4652
    FFV1_0( w_fp[449], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4653 OF 15495 ***
    // Wavefunction(s) for diagram number 4653
    // (none)
    // Amplitude(s) for diagram number 4653
    FFV1_0( w_fp[168], w_fp[473], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

    // *** DIAGRAM 4654 OF 15495 ***
    // Wavefunction(s) for diagram number 4654
    // (none)
    // Amplitude(s) for diagram number 4654
    FFV1_0( w_fp[256], w_fp[215], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4655 OF 15495 ***
    // Wavefunction(s) for diagram number 4655
    // (none)
    // Amplitude(s) for diagram number 4655
    VVV1_0( w_fp[488], w_fp[1], w_fp[231], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4656 OF 15495 ***
    // Wavefunction(s) for diagram number 4656
    // (none)
    // Amplitude(s) for diagram number 4656
    FFV1_0( w_fp[170], w_fp[473], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4657 OF 15495 ***
    // Wavefunction(s) for diagram number 4657
    // (none)
    // Amplitude(s) for diagram number 4657
    FFV1_0( w_fp[256], w_fp[227], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4658 OF 15495 ***
    // Wavefunction(s) for diagram number 4658
    // (none)
    // Amplitude(s) for diagram number 4658
    FFV1_0( w_fp[168], w_fp[215], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4659 OF 15495 ***
    // Wavefunction(s) for diagram number 4659
    // (none)
    // Amplitude(s) for diagram number 4659
    VVV1_0( w_fp[546], w_fp[232], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4660 OF 15495 ***
    // Wavefunction(s) for diagram number 4660
    // (none)
    // Amplitude(s) for diagram number 4660
    FFV1_0( w_fp[174], w_fp[226], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4661 OF 15495 ***
    // Wavefunction(s) for diagram number 4661
    // (none)
    // Amplitude(s) for diagram number 4661
    FFV1_0( w_fp[175], w_fp[215], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];

    // *** DIAGRAM 4662 OF 15495 ***
    // Wavefunction(s) for diagram number 4662
    // (none)
    // Amplitude(s) for diagram number 4662
    FFV1_0( w_fp[255], w_fp[543], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4663 OF 15495 ***
    // Wavefunction(s) for diagram number 4663
    // (none)
    // Amplitude(s) for diagram number 4663
    FFV1_0( w_fp[175], w_fp[543], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4664 OF 15495 ***
    // Wavefunction(s) for diagram number 4664
    // (none)
    // Amplitude(s) for diagram number 4664
    FFV1_0( w_fp[471], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4665 OF 15495 ***
    // Wavefunction(s) for diagram number 4665
    // (none)
    // Amplitude(s) for diagram number 4665
    FFV1_0( w_fp[471], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4666 OF 15495 ***
    // Wavefunction(s) for diagram number 4666
    // (none)
    // Amplitude(s) for diagram number 4666
    FFV1_0( w_fp[174], w_fp[473], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 4667 OF 15495 ***
    // Wavefunction(s) for diagram number 4667
    // (none)
    // Amplitude(s) for diagram number 4667
    FFV1_0( w_fp[255], w_fp[215], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4668 OF 15495 ***
    // Wavefunction(s) for diagram number 4668
    // (none)
    // Amplitude(s) for diagram number 4668
    VVV1_0( w_fp[451], w_fp[1], w_fp[232], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4669 OF 15495 ***
    // Wavefunction(s) for diagram number 4669
    // (none)
    // Amplitude(s) for diagram number 4669
    FFV1_0( w_fp[175], w_fp[473], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4670 OF 15495 ***
    // Wavefunction(s) for diagram number 4670
    // (none)
    // Amplitude(s) for diagram number 4670
    FFV1_0( w_fp[255], w_fp[226], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4671 OF 15495 ***
    // Wavefunction(s) for diagram number 4671
    // (none)
    // Amplitude(s) for diagram number 4671
    FFV1_0( w_fp[174], w_fp[215], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4672 OF 15495 ***
    // Wavefunction(s) for diagram number 4672
    // (none)
    // Amplitude(s) for diagram number 4672
    VVV1_0( w_fp[546], w_fp[228], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 4673 OF 15495 ***
    // Wavefunction(s) for diagram number 4673
    // (none)
    // Amplitude(s) for diagram number 4673
    FFV1_0( w_fp[3], w_fp[235], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4674 OF 15495 ***
    // Wavefunction(s) for diagram number 4674
    // (none)
    // Amplitude(s) for diagram number 4674
    FFV1_0( w_fp[184], w_fp[215], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4675 OF 15495 ***
    // Wavefunction(s) for diagram number 4675
    // (none)
    // Amplitude(s) for diagram number 4675
    FFV1_0( w_fp[3], w_fp[543], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4676 OF 15495 ***
    // Wavefunction(s) for diagram number 4676
    // (none)
    // Amplitude(s) for diagram number 4676
    FFV1_0( w_fp[184], w_fp[543], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4677 OF 15495 ***
    // Wavefunction(s) for diagram number 4677
    // (none)
    // Amplitude(s) for diagram number 4677
    FFV1_0( w_fp[479], w_fp[473], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];

    // *** DIAGRAM 4678 OF 15495 ***
    // Wavefunction(s) for diagram number 4678
    // (none)
    // Amplitude(s) for diagram number 4678
    FFV1_0( w_fp[479], w_fp[215], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4679 OF 15495 ***
    // Wavefunction(s) for diagram number 4679
    // (none)
    // Amplitude(s) for diagram number 4679
    FFV1_0( w_fp[479], w_fp[235], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 4680 OF 15495 ***
    // Wavefunction(s) for diagram number 4680
    // (none)
    // Amplitude(s) for diagram number 4680
    FFV1_0( w_fp[3], w_fp[473], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4681 OF 15495 ***
    // Wavefunction(s) for diagram number 4681
    // (none)
    // Amplitude(s) for diagram number 4681
    VVV1_0( w_fp[518], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 4682 OF 15495 ***
    // Wavefunction(s) for diagram number 4682
    // (none)
    // Amplitude(s) for diagram number 4682
    FFV1_0( w_fp[184], w_fp[473], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];

    // *** DIAGRAM 4683 OF 15495 ***
    // Wavefunction(s) for diagram number 4683
    // (none)
    // Amplitude(s) for diagram number 4683
    VVV1_0( w_fp[435], w_fp[359], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 4684 OF 15495 ***
    // Wavefunction(s) for diagram number 4684
    // (none)
    // Amplitude(s) for diagram number 4684
    FFV1_0( w_fp[3], w_fp[215], w_fp[588], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[590], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 4685 OF 15495 ***
    // Wavefunction(s) for diagram number 4685
    // (none)
    // Amplitude(s) for diagram number 4685
    VVVV1_0( w_fp[546], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVVV3_0( w_fp[546], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVVV4_0( w_fp[546], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4686 OF 15495 ***
    // Wavefunction(s) for diagram number 4686
    // (none)
    // Amplitude(s) for diagram number 4686
    VVV1_0( w_fp[154], w_fp[7], w_fp[561], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4687 OF 15495 ***
    // Wavefunction(s) for diagram number 4687
    // (none)
    // Amplitude(s) for diagram number 4687
    VVV1_0( w_fp[154], w_fp[6], w_fp[550], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4688 OF 15495 ***
    // Wavefunction(s) for diagram number 4688
    FFV1_1( w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[453] );
    // Amplitude(s) for diagram number 4688
    FFV1_0( w_fp[170], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];

    // *** DIAGRAM 4689 OF 15495 ***
    // Wavefunction(s) for diagram number 4689
    // (none)
    // Amplitude(s) for diagram number 4689
    FFV1_0( w_fp[170], w_fp[2], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4690 OF 15495 ***
    // Wavefunction(s) for diagram number 4690
    // (none)
    // Amplitude(s) for diagram number 4690
    FFV1_0( w_fp[171], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];

    // *** DIAGRAM 4691 OF 15495 ***
    // Wavefunction(s) for diagram number 4691
    // (none)
    // Amplitude(s) for diagram number 4691
    FFV1_0( w_fp[171], w_fp[2], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4692 OF 15495 ***
    // Wavefunction(s) for diagram number 4692
    // (none)
    // Amplitude(s) for diagram number 4692
    FFV1_0( w_fp[256], w_fp[531], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4693 OF 15495 ***
    // Wavefunction(s) for diagram number 4693
    // (none)
    // Amplitude(s) for diagram number 4693
    FFV1_0( w_fp[256], w_fp[532], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4694 OF 15495 ***
    // Wavefunction(s) for diagram number 4694
    FFV1_1( w_fp[530], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[527] );
    // Amplitude(s) for diagram number 4694
    FFV1_0( w_fp[170], w_fp[527], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4695 OF 15495 ***
    // Wavefunction(s) for diagram number 4695
    // (none)
    // Amplitude(s) for diagram number 4695
    FFV1_0( w_fp[170], w_fp[532], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4696 OF 15495 ***
    // Wavefunction(s) for diagram number 4696
    // (none)
    // Amplitude(s) for diagram number 4696
    FFV1_0( w_fp[171], w_fp[527], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4697 OF 15495 ***
    // Wavefunction(s) for diagram number 4697
    // (none)
    // Amplitude(s) for diagram number 4697
    FFV1_0( w_fp[171], w_fp[531], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4698 OF 15495 ***
    // Wavefunction(s) for diagram number 4698
    // (none)
    // Amplitude(s) for diagram number 4698
    FFV1_0( w_fp[256], w_fp[536], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];

    // *** DIAGRAM 4699 OF 15495 ***
    // Wavefunction(s) for diagram number 4699
    // (none)
    // Amplitude(s) for diagram number 4699
    FFV1_0( w_fp[256], w_fp[2], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4700 OF 15495 ***
    // Wavefunction(s) for diagram number 4700
    // (none)
    // Amplitude(s) for diagram number 4700
    VVVV1_0( w_fp[488], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4701 OF 15495 ***
    // Wavefunction(s) for diagram number 4701
    // (none)
    // Amplitude(s) for diagram number 4701
    VVV1_0( w_fp[154], w_fp[7], w_fp[571], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4702 OF 15495 ***
    // Wavefunction(s) for diagram number 4702
    // (none)
    // Amplitude(s) for diagram number 4702
    VVV1_0( w_fp[1], w_fp[154], w_fp[45], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4703 OF 15495 ***
    // Wavefunction(s) for diagram number 4703
    // (none)
    // Amplitude(s) for diagram number 4703
    FFV1_0( w_fp[171], w_fp[2], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4704 OF 15495 ***
    // Wavefunction(s) for diagram number 4704
    // (none)
    // Amplitude(s) for diagram number 4704
    FFV1_0( w_fp[171], w_fp[536], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 4705 OF 15495 ***
    // Wavefunction(s) for diagram number 4705
    // (none)
    // Amplitude(s) for diagram number 4705
    FFV1_0( w_fp[256], w_fp[537], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 4706 OF 15495 ***
    // Wavefunction(s) for diagram number 4706
    // (none)
    // Amplitude(s) for diagram number 4706
    FFV1_0( w_fp[256], w_fp[2], w_fp[509], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4707 OF 15495 ***
    // Wavefunction(s) for diagram number 4707
    // (none)
    // Amplitude(s) for diagram number 4707
    VVVV1_0( w_fp[437], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    VVVV3_0( w_fp[437], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    VVVV4_0( w_fp[437], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 4708 OF 15495 ***
    // Wavefunction(s) for diagram number 4708
    // (none)
    // Amplitude(s) for diagram number 4708
    VVV1_0( w_fp[154], w_fp[6], w_fp[575], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4709 OF 15495 ***
    // Wavefunction(s) for diagram number 4709
    // (none)
    // Amplitude(s) for diagram number 4709
    VVV1_0( w_fp[1], w_fp[154], w_fp[509], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 4710 OF 15495 ***
    // Wavefunction(s) for diagram number 4710
    // (none)
    // Amplitude(s) for diagram number 4710
    FFV1_0( w_fp[170], w_fp[2], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4711 OF 15495 ***
    // Wavefunction(s) for diagram number 4711
    // (none)
    // Amplitude(s) for diagram number 4711
    FFV1_0( w_fp[170], w_fp[537], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4712 OF 15495 ***
    // Wavefunction(s) for diagram number 4712
    // (none)
    // Amplitude(s) for diagram number 4712
    VVV1_0( w_fp[579], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    VVV1_0( w_fp[580], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0( w_fp[581], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4713 OF 15495 ***
    // Wavefunction(s) for diagram number 4713
    // (none)
    // Amplitude(s) for diagram number 4713
    FFV1_0( w_fp[171], w_fp[2], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4714 OF 15495 ***
    // Wavefunction(s) for diagram number 4714
    // (none)
    // Amplitude(s) for diagram number 4714
    VVV1_0( w_fp[582], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    VVV1_0( w_fp[583], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    VVV1_0( w_fp[584], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];

    // *** DIAGRAM 4715 OF 15495 ***
    // Wavefunction(s) for diagram number 4715
    // (none)
    // Amplitude(s) for diagram number 4715
    FFV1_0( w_fp[170], w_fp[2], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4716 OF 15495 ***
    // Wavefunction(s) for diagram number 4716
    // (none)
    // Amplitude(s) for diagram number 4716
    FFV1_0( w_fp[256], w_fp[2], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[445], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4717 OF 15495 ***
    // Wavefunction(s) for diagram number 4717
    // (none)
    // Amplitude(s) for diagram number 4717
    VVV1_0( w_fp[444], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0( w_fp[452], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    VVV1_0( w_fp[445], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4718 OF 15495 ***
    // Wavefunction(s) for diagram number 4718
    // (none)
    // Amplitude(s) for diagram number 4718
    VVV1_0( w_fp[546], w_fp[154], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4719 OF 15495 ***
    // Wavefunction(s) for diagram number 4719
    // (none)
    // Amplitude(s) for diagram number 4719
    FFV1_0( w_fp[168], w_fp[128], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4720 OF 15495 ***
    // Wavefunction(s) for diagram number 4720
    // (none)
    // Amplitude(s) for diagram number 4720
    FFV1_0( w_fp[112], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4721 OF 15495 ***
    // Wavefunction(s) for diagram number 4721
    // (none)
    // Amplitude(s) for diagram number 4721
    FFV1_0( w_fp[256], w_fp[530], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];

    // *** DIAGRAM 4722 OF 15495 ***
    // Wavefunction(s) for diagram number 4722
    // (none)
    // Amplitude(s) for diagram number 4722
    FFV1_0( w_fp[168], w_fp[530], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4723 OF 15495 ***
    // Wavefunction(s) for diagram number 4723
    // (none)
    // Amplitude(s) for diagram number 4723
    FFV1_0( w_fp[112], w_fp[530], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 4724 OF 15495 ***
    // Wavefunction(s) for diagram number 4724
    // (none)
    // Amplitude(s) for diagram number 4724
    FFV1_0( w_fp[449], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4725 OF 15495 ***
    // Wavefunction(s) for diagram number 4725
    // (none)
    // Amplitude(s) for diagram number 4725
    FFV1_0( w_fp[449], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];

    // *** DIAGRAM 4726 OF 15495 ***
    // Wavefunction(s) for diagram number 4726
    // (none)
    // Amplitude(s) for diagram number 4726
    FFV1_0( w_fp[256], w_fp[2], w_fp[472], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4727 OF 15495 ***
    // Wavefunction(s) for diagram number 4727
    // (none)
    // Amplitude(s) for diagram number 4727
    VVV1_0( w_fp[472], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4728 OF 15495 ***
    // Wavefunction(s) for diagram number 4728
    // (none)
    // Amplitude(s) for diagram number 4728
    FFV1_0( w_fp[256], w_fp[128], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4729 OF 15495 ***
    // Wavefunction(s) for diagram number 4729
    // (none)
    // Amplitude(s) for diagram number 4729
    VVV1_0( w_fp[435], w_fp[361], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 4730 OF 15495 ***
    // Wavefunction(s) for diagram number 4730
    // (none)
    // Amplitude(s) for diagram number 4730
    FFV1_0( w_fp[168], w_fp[2], w_fp[587], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[168], w_fp[2], w_fp[586], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[168], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4731 OF 15495 ***
    // Wavefunction(s) for diagram number 4731
    // (none)
    // Amplitude(s) for diagram number 4731
    VVVV1_0( w_fp[546], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV3_0( w_fp[546], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV4_0( w_fp[546], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4732 OF 15495 ***
    // Wavefunction(s) for diagram number 4732
    // (none)
    // Amplitude(s) for diagram number 4732
    VVV1_0( w_fp[150], w_fp[7], w_fp[549], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4733 OF 15495 ***
    // Wavefunction(s) for diagram number 4733
    // (none)
    // Amplitude(s) for diagram number 4733
    VVV1_0( w_fp[150], w_fp[5], w_fp[550], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4734 OF 15495 ***
    // Wavefunction(s) for diagram number 4734
    // (none)
    // Amplitude(s) for diagram number 4734
    FFV1_0( w_fp[175], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];

    // *** DIAGRAM 4735 OF 15495 ***
    // Wavefunction(s) for diagram number 4735
    // (none)
    // Amplitude(s) for diagram number 4735
    FFV1_0( w_fp[175], w_fp[2], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4736 OF 15495 ***
    // Wavefunction(s) for diagram number 4736
    // (none)
    // Amplitude(s) for diagram number 4736
    FFV1_0( w_fp[176], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];

    // *** DIAGRAM 4737 OF 15495 ***
    // Wavefunction(s) for diagram number 4737
    // (none)
    // Amplitude(s) for diagram number 4737
    FFV1_0( w_fp[176], w_fp[2], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4738 OF 15495 ***
    // Wavefunction(s) for diagram number 4738
    // (none)
    // Amplitude(s) for diagram number 4738
    FFV1_0( w_fp[255], w_fp[533], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4739 OF 15495 ***
    // Wavefunction(s) for diagram number 4739
    // (none)
    // Amplitude(s) for diagram number 4739
    FFV1_0( w_fp[255], w_fp[532], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4740 OF 15495 ***
    // Wavefunction(s) for diagram number 4740
    // (none)
    // Amplitude(s) for diagram number 4740
    FFV1_0( w_fp[175], w_fp[527], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4741 OF 15495 ***
    // Wavefunction(s) for diagram number 4741
    // (none)
    // Amplitude(s) for diagram number 4741
    FFV1_0( w_fp[175], w_fp[532], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4742 OF 15495 ***
    // Wavefunction(s) for diagram number 4742
    // (none)
    // Amplitude(s) for diagram number 4742
    FFV1_0( w_fp[176], w_fp[527], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4743 OF 15495 ***
    // Wavefunction(s) for diagram number 4743
    // (none)
    // Amplitude(s) for diagram number 4743
    FFV1_0( w_fp[176], w_fp[533], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4744 OF 15495 ***
    // Wavefunction(s) for diagram number 4744
    // (none)
    // Amplitude(s) for diagram number 4744
    FFV1_0( w_fp[255], w_fp[535], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];

    // *** DIAGRAM 4745 OF 15495 ***
    // Wavefunction(s) for diagram number 4745
    // (none)
    // Amplitude(s) for diagram number 4745
    FFV1_0( w_fp[255], w_fp[2], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4746 OF 15495 ***
    // Wavefunction(s) for diagram number 4746
    // (none)
    // Amplitude(s) for diagram number 4746
    VVVV1_0( w_fp[451], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV3_0( w_fp[451], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    VVVV4_0( w_fp[451], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4747 OF 15495 ***
    // Wavefunction(s) for diagram number 4747
    // (none)
    // Amplitude(s) for diagram number 4747
    VVV1_0( w_fp[150], w_fp[7], w_fp[553], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4748 OF 15495 ***
    // Wavefunction(s) for diagram number 4748
    // (none)
    // Amplitude(s) for diagram number 4748
    VVV1_0( w_fp[1], w_fp[150], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4749 OF 15495 ***
    // Wavefunction(s) for diagram number 4749
    // (none)
    // Amplitude(s) for diagram number 4749
    FFV1_0( w_fp[176], w_fp[2], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4750 OF 15495 ***
    // Wavefunction(s) for diagram number 4750
    // (none)
    // Amplitude(s) for diagram number 4750
    FFV1_0( w_fp[176], w_fp[535], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];

    // *** DIAGRAM 4751 OF 15495 ***
    // Wavefunction(s) for diagram number 4751
    // (none)
    // Amplitude(s) for diagram number 4751
    FFV1_0( w_fp[255], w_fp[537], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];

    // *** DIAGRAM 4752 OF 15495 ***
    // Wavefunction(s) for diagram number 4752
    // (none)
    // Amplitude(s) for diagram number 4752
    FFV1_0( w_fp[255], w_fp[2], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4753 OF 15495 ***
    // Wavefunction(s) for diagram number 4753
    // (none)
    // Amplitude(s) for diagram number 4753
    VVVV1_0( w_fp[437], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    VVVV3_0( w_fp[437], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    VVVV4_0( w_fp[437], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 4754 OF 15495 ***
    // Wavefunction(s) for diagram number 4754
    // (none)
    // Amplitude(s) for diagram number 4754
    VVV1_0( w_fp[150], w_fp[5], w_fp[575], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4755 OF 15495 ***
    // Wavefunction(s) for diagram number 4755
    // (none)
    // Amplitude(s) for diagram number 4755
    VVV1_0( w_fp[1], w_fp[150], w_fp[480], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 4756 OF 15495 ***
    // Wavefunction(s) for diagram number 4756
    // (none)
    // Amplitude(s) for diagram number 4756
    FFV1_0( w_fp[175], w_fp[2], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4757 OF 15495 ***
    // Wavefunction(s) for diagram number 4757
    // (none)
    // Amplitude(s) for diagram number 4757
    FFV1_0( w_fp[175], w_fp[537], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4758 OF 15495 ***
    // Wavefunction(s) for diagram number 4758
    // (none)
    // Amplitude(s) for diagram number 4758
    VVV1_0( w_fp[576], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    VVV1_0( w_fp[577], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[578], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4759 OF 15495 ***
    // Wavefunction(s) for diagram number 4759
    // (none)
    // Amplitude(s) for diagram number 4759
    FFV1_0( w_fp[176], w_fp[2], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4760 OF 15495 ***
    // Wavefunction(s) for diagram number 4760
    // (none)
    // Amplitude(s) for diagram number 4760
    VVV1_0( w_fp[582], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    VVV1_0( w_fp[583], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    VVV1_0( w_fp[584], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];

    // *** DIAGRAM 4761 OF 15495 ***
    // Wavefunction(s) for diagram number 4761
    // (none)
    // Amplitude(s) for diagram number 4761
    FFV1_0( w_fp[175], w_fp[2], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4762 OF 15495 ***
    // Wavefunction(s) for diagram number 4762
    // (none)
    // Amplitude(s) for diagram number 4762
    FFV1_0( w_fp[255], w_fp[2], w_fp[481], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[448], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4763 OF 15495 ***
    // Wavefunction(s) for diagram number 4763
    // (none)
    // Amplitude(s) for diagram number 4763
    VVV1_0( w_fp[481], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[510], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    VVV1_0( w_fp[448], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4764 OF 15495 ***
    // Wavefunction(s) for diagram number 4764
    // (none)
    // Amplitude(s) for diagram number 4764
    VVV1_0( w_fp[546], w_fp[150], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4765 OF 15495 ***
    // Wavefunction(s) for diagram number 4765
    // (none)
    // Amplitude(s) for diagram number 4765
    FFV1_0( w_fp[174], w_fp[122], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4766 OF 15495 ***
    // Wavefunction(s) for diagram number 4766
    // (none)
    // Amplitude(s) for diagram number 4766
    FFV1_0( w_fp[109], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4767 OF 15495 ***
    // Wavefunction(s) for diagram number 4767
    // (none)
    // Amplitude(s) for diagram number 4767
    FFV1_0( w_fp[255], w_fp[530], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];

    // *** DIAGRAM 4768 OF 15495 ***
    // Wavefunction(s) for diagram number 4768
    // (none)
    // Amplitude(s) for diagram number 4768
    FFV1_0( w_fp[174], w_fp[530], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4769 OF 15495 ***
    // Wavefunction(s) for diagram number 4769
    // (none)
    // Amplitude(s) for diagram number 4769
    FFV1_0( w_fp[109], w_fp[530], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];

    // *** DIAGRAM 4770 OF 15495 ***
    // Wavefunction(s) for diagram number 4770
    // (none)
    // Amplitude(s) for diagram number 4770
    FFV1_0( w_fp[471], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4771 OF 15495 ***
    // Wavefunction(s) for diagram number 4771
    // (none)
    // Amplitude(s) for diagram number 4771
    FFV1_0( w_fp[471], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];

    // *** DIAGRAM 4772 OF 15495 ***
    // Wavefunction(s) for diagram number 4772
    // (none)
    // Amplitude(s) for diagram number 4772
    FFV1_0( w_fp[255], w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4773 OF 15495 ***
    // Wavefunction(s) for diagram number 4773
    // (none)
    // Amplitude(s) for diagram number 4773
    VVV1_0( w_fp[474], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4774 OF 15495 ***
    // Wavefunction(s) for diagram number 4774
    // (none)
    // Amplitude(s) for diagram number 4774
    FFV1_0( w_fp[255], w_fp[122], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4775 OF 15495 ***
    // Wavefunction(s) for diagram number 4775
    // (none)
    // Amplitude(s) for diagram number 4775
    VVV1_0( w_fp[435], w_fp[360], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 4776 OF 15495 ***
    // Wavefunction(s) for diagram number 4776
    // (none)
    // Amplitude(s) for diagram number 4776
    FFV1_0( w_fp[174], w_fp[2], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[174], w_fp[2], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[174], w_fp[2], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4777 OF 15495 ***
    // Wavefunction(s) for diagram number 4777
    // (none)
    // Amplitude(s) for diagram number 4777
    VVVV1_0( w_fp[546], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV3_0( w_fp[546], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV4_0( w_fp[546], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 4778 OF 15495 ***
    // Wavefunction(s) for diagram number 4778
    // (none)
    // Amplitude(s) for diagram number 4778
    VVV1_0( w_fp[144], w_fp[6], w_fp[549], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 4779 OF 15495 ***
    // Wavefunction(s) for diagram number 4779
    // (none)
    // Amplitude(s) for diagram number 4779
    VVV1_0( w_fp[144], w_fp[5], w_fp[561], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 4780 OF 15495 ***
    // Wavefunction(s) for diagram number 4780
    // (none)
    // Amplitude(s) for diagram number 4780
    FFV1_0( w_fp[180], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];

    // *** DIAGRAM 4781 OF 15495 ***
    // Wavefunction(s) for diagram number 4781
    // (none)
    // Amplitude(s) for diagram number 4781
    FFV1_0( w_fp[180], w_fp[2], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4782 OF 15495 ***
    // Wavefunction(s) for diagram number 4782
    // (none)
    // Amplitude(s) for diagram number 4782
    FFV1_0( w_fp[181], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];

    // *** DIAGRAM 4783 OF 15495 ***
    // Wavefunction(s) for diagram number 4783
    // (none)
    // Amplitude(s) for diagram number 4783
    FFV1_0( w_fp[181], w_fp[2], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4784 OF 15495 ***
    // Wavefunction(s) for diagram number 4784
    // (none)
    // Amplitude(s) for diagram number 4784
    FFV1_0( w_fp[252], w_fp[533], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4785 OF 15495 ***
    // Wavefunction(s) for diagram number 4785
    // (none)
    // Amplitude(s) for diagram number 4785
    FFV1_0( w_fp[252], w_fp[531], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4786 OF 15495 ***
    // Wavefunction(s) for diagram number 4786
    // (none)
    // Amplitude(s) for diagram number 4786
    FFV1_0( w_fp[180], w_fp[527], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4787 OF 15495 ***
    // Wavefunction(s) for diagram number 4787
    // (none)
    // Amplitude(s) for diagram number 4787
    FFV1_0( w_fp[180], w_fp[531], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4788 OF 15495 ***
    // Wavefunction(s) for diagram number 4788
    // (none)
    // Amplitude(s) for diagram number 4788
    FFV1_0( w_fp[181], w_fp[527], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4789 OF 15495 ***
    // Wavefunction(s) for diagram number 4789
    // (none)
    // Amplitude(s) for diagram number 4789
    FFV1_0( w_fp[181], w_fp[533], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4790 OF 15495 ***
    // Wavefunction(s) for diagram number 4790
    // (none)
    // Amplitude(s) for diagram number 4790
    FFV1_0( w_fp[252], w_fp[535], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];

    // *** DIAGRAM 4791 OF 15495 ***
    // Wavefunction(s) for diagram number 4791
    // (none)
    // Amplitude(s) for diagram number 4791
    FFV1_0( w_fp[252], w_fp[2], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4792 OF 15495 ***
    // Wavefunction(s) for diagram number 4792
    // (none)
    // Amplitude(s) for diagram number 4792
    VVVV1_0( w_fp[451], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV3_0( w_fp[451], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVVV4_0( w_fp[451], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 4793 OF 15495 ***
    // Wavefunction(s) for diagram number 4793
    // (none)
    // Amplitude(s) for diagram number 4793
    VVV1_0( w_fp[144], w_fp[6], w_fp[553], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 4794 OF 15495 ***
    // Wavefunction(s) for diagram number 4794
    // (none)
    // Amplitude(s) for diagram number 4794
    VVV1_0( w_fp[1], w_fp[144], w_fp[514], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 4795 OF 15495 ***
    // Wavefunction(s) for diagram number 4795
    // (none)
    // Amplitude(s) for diagram number 4795
    FFV1_0( w_fp[181], w_fp[2], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4796 OF 15495 ***
    // Wavefunction(s) for diagram number 4796
    // (none)
    // Amplitude(s) for diagram number 4796
    FFV1_0( w_fp[181], w_fp[535], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];

    // *** DIAGRAM 4797 OF 15495 ***
    // Wavefunction(s) for diagram number 4797
    // (none)
    // Amplitude(s) for diagram number 4797
    FFV1_0( w_fp[252], w_fp[536], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];

    // *** DIAGRAM 4798 OF 15495 ***
    // Wavefunction(s) for diagram number 4798
    // (none)
    // Amplitude(s) for diagram number 4798
    FFV1_0( w_fp[252], w_fp[2], w_fp[497], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4799 OF 15495 ***
    // Wavefunction(s) for diagram number 4799
    // (none)
    // Amplitude(s) for diagram number 4799
    VVVV1_0( w_fp[488], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 4800 OF 15495 ***
    // Wavefunction(s) for diagram number 4800
    // (none)
    // Amplitude(s) for diagram number 4800
    VVV1_0( w_fp[144], w_fp[5], w_fp[571], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 453 );
    storeWf( wfs, w_cx, nevt, 527 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
