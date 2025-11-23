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
  diagramgroup67( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 22 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 106 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 117 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 149 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 173 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 189 );
    retrieveWf( wfs, w_cx, nevt, 192 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 230 );
    retrieveWf( wfs, w_cx, nevt, 231 );
    retrieveWf( wfs, w_cx, nevt, 232 );
    retrieveWf( wfs, w_cx, nevt, 233 );
    retrieveWf( wfs, w_cx, nevt, 234 );
    retrieveWf( wfs, w_cx, nevt, 235 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 247 );
    retrieveWf( wfs, w_cx, nevt, 250 );
    retrieveWf( wfs, w_cx, nevt, 251 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 262 );
    retrieveWf( wfs, w_cx, nevt, 280 );
    retrieveWf( wfs, w_cx, nevt, 353 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 366 );
    retrieveWf( wfs, w_cx, nevt, 431 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 491 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 590 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 597 );
    retrieveWf( wfs, w_cx, nevt, 599 );
    retrieveWf( wfs, w_cx, nevt, 601 );
    retrieveWf( wfs, w_cx, nevt, 603 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 605 );
    retrieveWf( wfs, w_cx, nevt, 607 );
    retrieveWf( wfs, w_cx, nevt, 608 );
    retrieveWf( wfs, w_cx, nevt, 610 );
    retrieveWf( wfs, w_cx, nevt, 620 );
    retrieveWf( wfs, w_cx, nevt, 621 );
    retrieveWf( wfs, w_cx, nevt, 623 );
    retrieveWf( wfs, w_cx, nevt, 630 );
    retrieveWf( wfs, w_cx, nevt, 631 );
    retrieveWf( wfs, w_cx, nevt, 632 );
    retrieveWf( wfs, w_cx, nevt, 633 );
    retrieveWf( wfs, w_cx, nevt, 634 );
    retrieveWf( wfs, w_cx, nevt, 635 );
    retrieveWf( wfs, w_cx, nevt, 658 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 674 );
    retrieveWf( wfs, w_cx, nevt, 690 );
    retrieveWf( wfs, w_cx, nevt, 691 );
    retrieveWf( wfs, w_cx, nevt, 692 );
    retrieveWf( wfs, w_cx, nevt, 693 );
    retrieveWf( wfs, w_cx, nevt, 694 );
    retrieveWf( wfs, w_cx, nevt, 695 );
    retrieveWf( wfs, w_cx, nevt, 746 );
#endif
#endif

    // *** DIAGRAM 13201 OF 15495 ***
    // Wavefunction(s) for diagram number 13201
    // (none)
    // Amplitude(s) for diagram number 13201
    VVV1_0( w_fp[591], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[539], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    VVV1_0( w_fp[531], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 13202 OF 15495 ***
    // Wavefunction(s) for diagram number 13202
    // (none)
    // Amplitude(s) for diagram number 13202
    FFV1_0( w_fp[44], w_fp[476], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    FFV1_0( w_fp[192], w_fp[476], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    FFV1_0( w_fp[75], w_fp[476], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];

    // *** DIAGRAM 13203 OF 15495 ***
    // Wavefunction(s) for diagram number 13203
    // (none)
    // Amplitude(s) for diagram number 13203
    VVV1_0( w_fp[0], w_fp[108], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[0], w_fp[431], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[0], w_fp[280], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13204 OF 15495 ***
    // Wavefunction(s) for diagram number 13204
    // (none)
    // Amplitude(s) for diagram number 13204
    FFV1_0( w_fp[3], w_fp[197], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[162], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[635], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[634], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[633], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[632], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[631], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[630], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 13205 OF 15495 ***
    // Wavefunction(s) for diagram number 13205
    // (none)
    // Amplitude(s) for diagram number 13205
    FFV1_0( w_fp[355], w_fp[690], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] -= amp_sv[0];

    // *** DIAGRAM 13206 OF 15495 ***
    // Wavefunction(s) for diagram number 13206
    // (none)
    // Amplitude(s) for diagram number 13206
    FFV1_0( w_fp[355], w_fp[691], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] -= amp_sv[0];

    // *** DIAGRAM 13207 OF 15495 ***
    // Wavefunction(s) for diagram number 13207
    FFV1_1( w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[75] );
    // Amplitude(s) for diagram number 13207
    FFV1_0( w_fp[216], w_fp[75], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] -= amp_sv[0];

    // *** DIAGRAM 13208 OF 15495 ***
    // Wavefunction(s) for diagram number 13208
    // (none)
    // Amplitude(s) for diagram number 13208
    FFV1_0( w_fp[216], w_fp[691], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] -= amp_sv[0];

    // *** DIAGRAM 13209 OF 15495 ***
    // Wavefunction(s) for diagram number 13209
    // (none)
    // Amplitude(s) for diagram number 13209
    FFV1_0( w_fp[198], w_fp[75], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] -= amp_sv[0];

    // *** DIAGRAM 13210 OF 15495 ***
    // Wavefunction(s) for diagram number 13210
    // (none)
    // Amplitude(s) for diagram number 13210
    FFV1_0( w_fp[198], w_fp[690], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[613] -= amp_sv[0];

    // *** DIAGRAM 13211 OF 15495 ***
    // Wavefunction(s) for diagram number 13211
    // (none)
    // Amplitude(s) for diagram number 13211
    FFV1_0( w_fp[586], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[646] -= amp_sv[0];

    // *** DIAGRAM 13212 OF 15495 ***
    // Wavefunction(s) for diagram number 13212
    // (none)
    // Amplitude(s) for diagram number 13212
    FFV1_0( w_fp[560], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[640] -= amp_sv[0];

    // *** DIAGRAM 13213 OF 15495 ***
    // Wavefunction(s) for diagram number 13213
    // (none)
    // Amplitude(s) for diagram number 13213
    FFV1_0( w_fp[160], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[692] -= amp_sv[0];

    // *** DIAGRAM 13214 OF 15495 ***
    // Wavefunction(s) for diagram number 13214
    // (none)
    // Amplitude(s) for diagram number 13214
    FFV1_0( w_fp[560], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[682] -= amp_sv[0];

    // *** DIAGRAM 13215 OF 15495 ***
    // Wavefunction(s) for diagram number 13215
    // (none)
    // Amplitude(s) for diagram number 13215
    FFV1_0( w_fp[160], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 13216 OF 15495 ***
    // Wavefunction(s) for diagram number 13216
    // (none)
    // Amplitude(s) for diagram number 13216
    FFV1_0( w_fp[586], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 13217 OF 15495 ***
    // Wavefunction(s) for diagram number 13217
    FFV1_1( w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[192] );
    // Amplitude(s) for diagram number 13217
    FFV1_0( w_fp[216], w_fp[192], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 13218 OF 15495 ***
    // Wavefunction(s) for diagram number 13218
    // (none)
    // Amplitude(s) for diagram number 13218
    FFV1_0( w_fp[588], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[643] -= amp_sv[0];

    // *** DIAGRAM 13219 OF 15495 ***
    // Wavefunction(s) for diagram number 13219
    // (none)
    // Amplitude(s) for diagram number 13219
    FFV1_0( w_fp[198], w_fp[192], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[627] -= amp_sv[0];

    // *** DIAGRAM 13220 OF 15495 ***
    // Wavefunction(s) for diagram number 13220
    // (none)
    // Amplitude(s) for diagram number 13220
    FFV1_0( w_fp[541], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[637] -= amp_sv[0];

    // *** DIAGRAM 13221 OF 15495 ***
    // Wavefunction(s) for diagram number 13221
    // (none)
    // Amplitude(s) for diagram number 13221
    FFV1_0( w_fp[674], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] -= amp_sv[0];

    // *** DIAGRAM 13222 OF 15495 ***
    // Wavefunction(s) for diagram number 13222
    // (none)
    // Amplitude(s) for diagram number 13222
    FFV1_0( w_fp[355], w_fp[693], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] -= amp_sv[0];

    // *** DIAGRAM 13223 OF 15495 ***
    // Wavefunction(s) for diagram number 13223
    // (none)
    // Amplitude(s) for diagram number 13223
    FFV1_0( w_fp[674], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13224 OF 15495 ***
    // Wavefunction(s) for diagram number 13224
    // (none)
    // Amplitude(s) for diagram number 13224
    FFV1_0( w_fp[355], w_fp[694], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 13225 OF 15495 ***
    // Wavefunction(s) for diagram number 13225
    // (none)
    // Amplitude(s) for diagram number 13225
    FFV1_0( w_fp[198], w_fp[693], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[673] -= amp_sv[0];

    // *** DIAGRAM 13226 OF 15495 ***
    // Wavefunction(s) for diagram number 13226
    // (none)
    // Amplitude(s) for diagram number 13226
    FFV1_0( w_fp[541], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 13227 OF 15495 ***
    // Wavefunction(s) for diagram number 13227
    // (none)
    // Amplitude(s) for diagram number 13227
    FFV1_0( w_fp[216], w_fp[694], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[697] -= amp_sv[0];

    // *** DIAGRAM 13228 OF 15495 ***
    // Wavefunction(s) for diagram number 13228
    // (none)
    // Amplitude(s) for diagram number 13228
    FFV1_0( w_fp[588], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 13229 OF 15495 ***
    // Wavefunction(s) for diagram number 13229
    // (none)
    // Amplitude(s) for diagram number 13229
    FFV1_0( w_fp[355], w_fp[542], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13230 OF 15495 ***
    // Wavefunction(s) for diagram number 13230
    // (none)
    // Amplitude(s) for diagram number 13230
    FFV1_0( w_fp[196], w_fp[542], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 13231 OF 15495 ***
    // Wavefunction(s) for diagram number 13231
    // (none)
    // Amplitude(s) for diagram number 13231
    FFV1_0( w_fp[243], w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13232 OF 15495 ***
    // Wavefunction(s) for diagram number 13232
    // (none)
    // Amplitude(s) for diagram number 13232
    FFV1_0( w_fp[536], w_fp[473], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13233 OF 15495 ***
    // Wavefunction(s) for diagram number 13233
    // (none)
    // Amplitude(s) for diagram number 13233
    FFV1_0( w_fp[536], w_fp[215], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13234 OF 15495 ***
    // Wavefunction(s) for diagram number 13234
    // (none)
    // Amplitude(s) for diagram number 13234
    FFV1_0( w_fp[536], w_fp[235], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13235 OF 15495 ***
    // Wavefunction(s) for diagram number 13235
    // (none)
    // Amplitude(s) for diagram number 13235
    FFV1_0( w_fp[196], w_fp[473], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];

    // *** DIAGRAM 13236 OF 15495 ***
    // Wavefunction(s) for diagram number 13236
    // (none)
    // Amplitude(s) for diagram number 13236
    FFV1_0( w_fp[355], w_fp[215], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13237 OF 15495 ***
    // Wavefunction(s) for diagram number 13237
    // (none)
    // Amplitude(s) for diagram number 13237
    VVV1_0( w_fp[559], w_fp[1], w_fp[230], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13238 OF 15495 ***
    // Wavefunction(s) for diagram number 13238
    // (none)
    // Amplitude(s) for diagram number 13238
    FFV1_0( w_fp[243], w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13239 OF 15495 ***
    // Wavefunction(s) for diagram number 13239
    // (none)
    // Amplitude(s) for diagram number 13239
    FFV1_0( w_fp[355], w_fp[235], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13240 OF 15495 ***
    // Wavefunction(s) for diagram number 13240
    // (none)
    // Amplitude(s) for diagram number 13240
    VVV1_0( w_fp[0], w_fp[359], w_fp[230], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13241 OF 15495 ***
    // Wavefunction(s) for diagram number 13241
    // (none)
    // Amplitude(s) for diagram number 13241
    FFV1_0( w_fp[196], w_fp[215], w_fp[106], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[215], w_fp[117], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[215], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13242 OF 15495 ***
    // Wavefunction(s) for diagram number 13242
    // (none)
    // Amplitude(s) for diagram number 13242
    FFV1_0( w_fp[256], w_fp[692], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= amp_sv[0];

    // *** DIAGRAM 13243 OF 15495 ***
    // Wavefunction(s) for diagram number 13243
    // (none)
    // Amplitude(s) for diagram number 13243
    FFV1_0( w_fp[256], w_fp[691], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] -= amp_sv[0];

    // *** DIAGRAM 13244 OF 15495 ***
    // Wavefunction(s) for diagram number 13244
    // (none)
    // Amplitude(s) for diagram number 13244
    FFV1_0( w_fp[218], w_fp[75], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= amp_sv[0];

    // *** DIAGRAM 13245 OF 15495 ***
    // Wavefunction(s) for diagram number 13245
    // (none)
    // Amplitude(s) for diagram number 13245
    FFV1_0( w_fp[218], w_fp[691], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];

    // *** DIAGRAM 13246 OF 15495 ***
    // Wavefunction(s) for diagram number 13246
    // (none)
    // Amplitude(s) for diagram number 13246
    FFV1_0( w_fp[170], w_fp[75], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= amp_sv[0];

    // *** DIAGRAM 13247 OF 15495 ***
    // Wavefunction(s) for diagram number 13247
    // (none)
    // Amplitude(s) for diagram number 13247
    FFV1_0( w_fp[170], w_fp[692], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] -= amp_sv[0];

    // *** DIAGRAM 13248 OF 15495 ***
    // Wavefunction(s) for diagram number 13248
    // (none)
    // Amplitude(s) for diagram number 13248
    FFV1_0( w_fp[509], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[644] -= amp_sv[0];

    // *** DIAGRAM 13249 OF 15495 ***
    // Wavefunction(s) for diagram number 13249
    // (none)
    // Amplitude(s) for diagram number 13249
    FFV1_0( w_fp[538], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] -= amp_sv[0];

    // *** DIAGRAM 13250 OF 15495 ***
    // Wavefunction(s) for diagram number 13250
    // (none)
    // Amplitude(s) for diagram number 13250
    FFV1_0( w_fp[660], w_fp[225], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[668] -= amp_sv[0];

    // *** DIAGRAM 13251 OF 15495 ***
    // Wavefunction(s) for diagram number 13251
    // (none)
    // Amplitude(s) for diagram number 13251
    FFV1_0( w_fp[538], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13252 OF 15495 ***
    // Wavefunction(s) for diagram number 13252
    // (none)
    // Amplitude(s) for diagram number 13252
    FFV1_0( w_fp[660], w_fp[227], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 13253 OF 15495 ***
    // Wavefunction(s) for diagram number 13253
    // (none)
    // Amplitude(s) for diagram number 13253
    FFV1_0( w_fp[509], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 13254 OF 15495 ***
    // Wavefunction(s) for diagram number 13254
    // (none)
    // Amplitude(s) for diagram number 13254
    FFV1_0( w_fp[218], w_fp[192], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 13255 OF 15495 ***
    // Wavefunction(s) for diagram number 13255
    // (none)
    // Amplitude(s) for diagram number 13255
    FFV1_0( w_fp[481], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] -= amp_sv[0];

    // *** DIAGRAM 13256 OF 15495 ***
    // Wavefunction(s) for diagram number 13256
    // (none)
    // Amplitude(s) for diagram number 13256
    FFV1_0( w_fp[170], w_fp[192], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] -= amp_sv[0];

    // *** DIAGRAM 13257 OF 15495 ***
    // Wavefunction(s) for diagram number 13257
    // (none)
    // Amplitude(s) for diagram number 13257
    FFV1_0( w_fp[535], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[631] -= amp_sv[0];

    // *** DIAGRAM 13258 OF 15495 ***
    // Wavefunction(s) for diagram number 13258
    // (none)
    // Amplitude(s) for diagram number 13258
    FFV1_0( w_fp[658], w_fp[225], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[666] -= amp_sv[0];

    // *** DIAGRAM 13259 OF 15495 ***
    // Wavefunction(s) for diagram number 13259
    // (none)
    // Amplitude(s) for diagram number 13259
    FFV1_0( w_fp[256], w_fp[695], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 13260 OF 15495 ***
    // Wavefunction(s) for diagram number 13260
    // (none)
    // Amplitude(s) for diagram number 13260
    FFV1_0( w_fp[658], w_fp[227], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13261 OF 15495 ***
    // Wavefunction(s) for diagram number 13261
    // (none)
    // Amplitude(s) for diagram number 13261
    FFV1_0( w_fp[256], w_fp[694], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[698] -= amp_sv[0];

    // *** DIAGRAM 13262 OF 15495 ***
    // Wavefunction(s) for diagram number 13262
    // (none)
    // Amplitude(s) for diagram number 13262
    FFV1_0( w_fp[170], w_fp[695], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[649] -= amp_sv[0];

    // *** DIAGRAM 13263 OF 15495 ***
    // Wavefunction(s) for diagram number 13263
    // (none)
    // Amplitude(s) for diagram number 13263
    FFV1_0( w_fp[535], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 13264 OF 15495 ***
    // Wavefunction(s) for diagram number 13264
    // (none)
    // Amplitude(s) for diagram number 13264
    FFV1_0( w_fp[218], w_fp[694], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] -= amp_sv[0];

    // *** DIAGRAM 13265 OF 15495 ***
    // Wavefunction(s) for diagram number 13265
    // (none)
    // Amplitude(s) for diagram number 13265
    FFV1_0( w_fp[481], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 13266 OF 15495 ***
    // Wavefunction(s) for diagram number 13266
    // (none)
    // Amplitude(s) for diagram number 13266
    FFV1_0( w_fp[256], w_fp[542], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13267 OF 15495 ***
    // Wavefunction(s) for diagram number 13267
    // (none)
    // Amplitude(s) for diagram number 13267
    FFV1_0( w_fp[168], w_fp[542], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 13268 OF 15495 ***
    // Wavefunction(s) for diagram number 13268
    // (none)
    // Amplitude(s) for diagram number 13268
    FFV1_0( w_fp[200], w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13269 OF 15495 ***
    // Wavefunction(s) for diagram number 13269
    // (none)
    // Amplitude(s) for diagram number 13269
    FFV1_0( w_fp[527], w_fp[473], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13270 OF 15495 ***
    // Wavefunction(s) for diagram number 13270
    // (none)
    // Amplitude(s) for diagram number 13270
    FFV1_0( w_fp[527], w_fp[215], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13271 OF 15495 ***
    // Wavefunction(s) for diagram number 13271
    // (none)
    // Amplitude(s) for diagram number 13271
    FFV1_0( w_fp[527], w_fp[234], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13272 OF 15495 ***
    // Wavefunction(s) for diagram number 13272
    // (none)
    // Amplitude(s) for diagram number 13272
    FFV1_0( w_fp[168], w_fp[473], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

    // *** DIAGRAM 13273 OF 15495 ***
    // Wavefunction(s) for diagram number 13273
    // (none)
    // Amplitude(s) for diagram number 13273
    FFV1_0( w_fp[256], w_fp[215], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13274 OF 15495 ***
    // Wavefunction(s) for diagram number 13274
    // (none)
    // Amplitude(s) for diagram number 13274
    VVV1_0( w_fp[576], w_fp[1], w_fp[231], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13275 OF 15495 ***
    // Wavefunction(s) for diagram number 13275
    // (none)
    // Amplitude(s) for diagram number 13275
    FFV1_0( w_fp[200], w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13276 OF 15495 ***
    // Wavefunction(s) for diagram number 13276
    // (none)
    // Amplitude(s) for diagram number 13276
    FFV1_0( w_fp[256], w_fp[234], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13277 OF 15495 ***
    // Wavefunction(s) for diagram number 13277
    // (none)
    // Amplitude(s) for diagram number 13277
    VVV1_0( w_fp[0], w_fp[362], w_fp[231], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13278 OF 15495 ***
    // Wavefunction(s) for diagram number 13278
    // (none)
    // Amplitude(s) for diagram number 13278
    FFV1_0( w_fp[168], w_fp[215], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13279 OF 15495 ***
    // Wavefunction(s) for diagram number 13279
    // (none)
    // Amplitude(s) for diagram number 13279
    FFV1_0( w_fp[255], w_fp[692], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];

    // *** DIAGRAM 13280 OF 15495 ***
    // Wavefunction(s) for diagram number 13280
    // (none)
    // Amplitude(s) for diagram number 13280
    FFV1_0( w_fp[255], w_fp[690], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] -= amp_sv[0];

    // *** DIAGRAM 13281 OF 15495 ***
    // Wavefunction(s) for diagram number 13281
    // (none)
    // Amplitude(s) for diagram number 13281
    FFV1_0( w_fp[202], w_fp[75], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[602] -= amp_sv[0];

    // *** DIAGRAM 13282 OF 15495 ***
    // Wavefunction(s) for diagram number 13282
    // (none)
    // Amplitude(s) for diagram number 13282
    FFV1_0( w_fp[202], w_fp[690], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= amp_sv[0];

    // *** DIAGRAM 13283 OF 15495 ***
    // Wavefunction(s) for diagram number 13283
    // (none)
    // Amplitude(s) for diagram number 13283
    FFV1_0( w_fp[175], w_fp[75], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];

    // *** DIAGRAM 13284 OF 15495 ***
    // Wavefunction(s) for diagram number 13284
    // (none)
    // Amplitude(s) for diagram number 13284
    FFV1_0( w_fp[175], w_fp[692], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] -= amp_sv[0];

    // *** DIAGRAM 13285 OF 15495 ***
    // Wavefunction(s) for diagram number 13285
    // (none)
    // Amplitude(s) for diagram number 13285
    FFV1_0( w_fp[595], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[638] -= amp_sv[0];

    // *** DIAGRAM 13286 OF 15495 ***
    // Wavefunction(s) for diagram number 13286
    // (none)
    // Amplitude(s) for diagram number 13286
    FFV1_0( w_fp[577], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] -= amp_sv[0];

    // *** DIAGRAM 13287 OF 15495 ***
    // Wavefunction(s) for diagram number 13287
    // (none)
    // Amplitude(s) for diagram number 13287
    FFV1_0( w_fp[173], w_fp[225], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[662] -= amp_sv[0];

    // *** DIAGRAM 13288 OF 15495 ***
    // Wavefunction(s) for diagram number 13288
    // (none)
    // Amplitude(s) for diagram number 13288
    FFV1_0( w_fp[577], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 13289 OF 15495 ***
    // Wavefunction(s) for diagram number 13289
    // (none)
    // Amplitude(s) for diagram number 13289
    FFV1_0( w_fp[173], w_fp[226], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 13290 OF 15495 ***
    // Wavefunction(s) for diagram number 13290
    // (none)
    // Amplitude(s) for diagram number 13290
    FFV1_0( w_fp[595], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 13291 OF 15495 ***
    // Wavefunction(s) for diagram number 13291
    // (none)
    // Amplitude(s) for diagram number 13291
    FFV1_0( w_fp[202], w_fp[192], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] -= amp_sv[0];

    // *** DIAGRAM 13292 OF 15495 ***
    // Wavefunction(s) for diagram number 13292
    // (none)
    // Amplitude(s) for diagram number 13292
    FFV1_0( w_fp[545], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] -= amp_sv[0];

    // *** DIAGRAM 13293 OF 15495 ***
    // Wavefunction(s) for diagram number 13293
    // (none)
    // Amplitude(s) for diagram number 13293
    FFV1_0( w_fp[175], w_fp[192], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= amp_sv[0];

    // *** DIAGRAM 13294 OF 15495 ***
    // Wavefunction(s) for diagram number 13294
    // (none)
    // Amplitude(s) for diagram number 13294
    FFV1_0( w_fp[530], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] -= amp_sv[0];

    // *** DIAGRAM 13295 OF 15495 ***
    // Wavefunction(s) for diagram number 13295
    // (none)
    // Amplitude(s) for diagram number 13295
    FFV1_0( w_fp[189], w_fp[225], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[660] -= amp_sv[0];

    // *** DIAGRAM 13296 OF 15495 ***
    // Wavefunction(s) for diagram number 13296
    // (none)
    // Amplitude(s) for diagram number 13296
    FFV1_0( w_fp[255], w_fp[695], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 13297 OF 15495 ***
    // Wavefunction(s) for diagram number 13297
    // (none)
    // Amplitude(s) for diagram number 13297
    FFV1_0( w_fp[189], w_fp[226], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13298 OF 15495 ***
    // Wavefunction(s) for diagram number 13298
    // (none)
    // Amplitude(s) for diagram number 13298
    FFV1_0( w_fp[255], w_fp[693], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 13299 OF 15495 ***
    // Wavefunction(s) for diagram number 13299
    // (none)
    // Amplitude(s) for diagram number 13299
    FFV1_0( w_fp[175], w_fp[695], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] -= amp_sv[0];

    // *** DIAGRAM 13300 OF 15495 ***
    // Wavefunction(s) for diagram number 13300
    // (none)
    // Amplitude(s) for diagram number 13300
    FFV1_0( w_fp[530], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 13301 OF 15495 ***
    // Wavefunction(s) for diagram number 13301
    // (none)
    // Amplitude(s) for diagram number 13301
    FFV1_0( w_fp[202], w_fp[693], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= amp_sv[0];

    // *** DIAGRAM 13302 OF 15495 ***
    // Wavefunction(s) for diagram number 13302
    // (none)
    // Amplitude(s) for diagram number 13302
    FFV1_0( w_fp[545], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 13303 OF 15495 ***
    // Wavefunction(s) for diagram number 13303
    // (none)
    // Amplitude(s) for diagram number 13303
    FFV1_0( w_fp[255], w_fp[542], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13304 OF 15495 ***
    // Wavefunction(s) for diagram number 13304
    // (none)
    // Amplitude(s) for diagram number 13304
    FFV1_0( w_fp[174], w_fp[542], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];

    // *** DIAGRAM 13305 OF 15495 ***
    // Wavefunction(s) for diagram number 13305
    // (none)
    // Amplitude(s) for diagram number 13305
    FFV1_0( w_fp[149], w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13306 OF 15495 ***
    // Wavefunction(s) for diagram number 13306
    // (none)
    // Amplitude(s) for diagram number 13306
    FFV1_0( w_fp[594], w_fp[473], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13307 OF 15495 ***
    // Wavefunction(s) for diagram number 13307
    // (none)
    // Amplitude(s) for diagram number 13307
    FFV1_0( w_fp[594], w_fp[215], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13308 OF 15495 ***
    // Wavefunction(s) for diagram number 13308
    // (none)
    // Amplitude(s) for diagram number 13308
    FFV1_0( w_fp[594], w_fp[233], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13309 OF 15495 ***
    // Wavefunction(s) for diagram number 13309
    // (none)
    // Amplitude(s) for diagram number 13309
    FFV1_0( w_fp[174], w_fp[473], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 13310 OF 15495 ***
    // Wavefunction(s) for diagram number 13310
    // (none)
    // Amplitude(s) for diagram number 13310
    FFV1_0( w_fp[255], w_fp[215], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13311 OF 15495 ***
    // Wavefunction(s) for diagram number 13311
    // (none)
    // Amplitude(s) for diagram number 13311
    VVV1_0( w_fp[514], w_fp[1], w_fp[232], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13312 OF 15495 ***
    // Wavefunction(s) for diagram number 13312
    // (none)
    // Amplitude(s) for diagram number 13312
    FFV1_0( w_fp[149], w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13313 OF 15495 ***
    // Wavefunction(s) for diagram number 13313
    // (none)
    // Amplitude(s) for diagram number 13313
    FFV1_0( w_fp[255], w_fp[233], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13314 OF 15495 ***
    // Wavefunction(s) for diagram number 13314
    // (none)
    // Amplitude(s) for diagram number 13314
    VVV1_0( w_fp[0], w_fp[254], w_fp[232], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13315 OF 15495 ***
    // Wavefunction(s) for diagram number 13315
    // (none)
    // Amplitude(s) for diagram number 13315
    FFV1_0( w_fp[174], w_fp[215], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13316 OF 15495 ***
    // Wavefunction(s) for diagram number 13316
    // (none)
    // Amplitude(s) for diagram number 13316
    VVV1_0( w_fp[254], w_fp[6], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13317 OF 15495 ***
    // Wavefunction(s) for diagram number 13317
    // (none)
    // Amplitude(s) for diagram number 13317
    FFV1_0( w_fp[3], w_fp[691], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];

    // *** DIAGRAM 13318 OF 15495 ***
    // Wavefunction(s) for diagram number 13318
    // (none)
    // Amplitude(s) for diagram number 13318
    FFV1_0( w_fp[221], w_fp[75], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13319 OF 15495 ***
    // Wavefunction(s) for diagram number 13319
    // (none)
    // Amplitude(s) for diagram number 13319
    FFV1_0( w_fp[221], w_fp[691], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13320 OF 15495 ***
    // Wavefunction(s) for diagram number 13320
    // (none)
    // Amplitude(s) for diagram number 13320
    FFV1_0( w_fp[3], w_fp[75], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 13321 OF 15495 ***
    // Wavefunction(s) for diagram number 13321
    // (none)
    // Amplitude(s) for diagram number 13321
    VVV1_0( w_fp[1], w_fp[68], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13322 OF 15495 ***
    // Wavefunction(s) for diagram number 13322
    // (none)
    // Amplitude(s) for diagram number 13322
    FFV1_0( w_fp[3], w_fp[542], w_fp[247], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[251], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[250], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13323 OF 15495 ***
    // Wavefunction(s) for diagram number 13323
    // (none)
    // Amplitude(s) for diagram number 13323
    FFV1_0( w_fp[599], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];

    // *** DIAGRAM 13324 OF 15495 ***
    // Wavefunction(s) for diagram number 13324
    // (none)
    // Amplitude(s) for diagram number 13324
    FFV1_0( w_fp[3], w_fp[473], w_fp[597], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13325 OF 15495 ***
    // Wavefunction(s) for diagram number 13325
    // (none)
    // Amplitude(s) for diagram number 13325
    VVVV1_0( w_fp[514], w_fp[1], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    VVVV3_0( w_fp[514], w_fp[1], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV4_0( w_fp[514], w_fp[1], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 13326 OF 15495 ***
    // Wavefunction(s) for diagram number 13326
    // (none)
    // Amplitude(s) for diagram number 13326
    VVV1_0( w_fp[228], w_fp[6], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 13327 OF 15495 ***
    // Wavefunction(s) for diagram number 13327
    // (none)
    // Amplitude(s) for diagram number 13327
    VVV1_0( w_fp[1], w_fp[228], w_fp[597], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 13328 OF 15495 ***
    // Wavefunction(s) for diagram number 13328
    // (none)
    // Amplitude(s) for diagram number 13328
    FFV1_0( w_fp[3], w_fp[227], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13329 OF 15495 ***
    // Wavefunction(s) for diagram number 13329
    // (none)
    // Amplitude(s) for diagram number 13329
    FFV1_0( w_fp[599], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 13330 OF 15495 ***
    // Wavefunction(s) for diagram number 13330
    // (none)
    // Amplitude(s) for diagram number 13330
    FFV1_0( w_fp[221], w_fp[192], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13331 OF 15495 ***
    // Wavefunction(s) for diagram number 13331
    // (none)
    // Amplitude(s) for diagram number 13331
    FFV1_0( w_fp[603], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13332 OF 15495 ***
    // Wavefunction(s) for diagram number 13332
    // (none)
    // Amplitude(s) for diagram number 13332
    FFV1_0( w_fp[3], w_fp[192], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 13333 OF 15495 ***
    // Wavefunction(s) for diagram number 13333
    // (none)
    // Amplitude(s) for diagram number 13333
    FFV1_0( w_fp[3], w_fp[473], w_fp[601], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13334 OF 15495 ***
    // Wavefunction(s) for diagram number 13334
    // (none)
    // Amplitude(s) for diagram number 13334
    VVVV1_0( w_fp[0], w_fp[254], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[254], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[254], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 13335 OF 15495 ***
    // Wavefunction(s) for diagram number 13335
    // (none)
    // Amplitude(s) for diagram number 13335
    VVV1_0( w_fp[228], w_fp[6], w_fp[182], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 13336 OF 15495 ***
    // Wavefunction(s) for diagram number 13336
    // (none)
    // Amplitude(s) for diagram number 13336
    VVV1_0( w_fp[254], w_fp[6], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 13337 OF 15495 ***
    // Wavefunction(s) for diagram number 13337
    // (none)
    // Amplitude(s) for diagram number 13337
    FFV1_0( w_fp[3], w_fp[227], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13338 OF 15495 ***
    // Wavefunction(s) for diagram number 13338
    // (none)
    // Amplitude(s) for diagram number 13338
    FFV1_0( w_fp[3], w_fp[694], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 13339 OF 15495 ***
    // Wavefunction(s) for diagram number 13339
    // (none)
    // Amplitude(s) for diagram number 13339
    VVVV1_0( w_fp[0], w_fp[1], w_fp[228], w_fp[68], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[228], w_fp[68], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[228], w_fp[68], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 13340 OF 15495 ***
    // Wavefunction(s) for diagram number 13340
    // (none)
    // Amplitude(s) for diagram number 13340
    VVV1_0( w_fp[1], w_fp[68], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 13341 OF 15495 ***
    // Wavefunction(s) for diagram number 13341
    // (none)
    // Amplitude(s) for diagram number 13341
    VVV1_0( w_fp[1], w_fp[228], w_fp[601], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 13342 OF 15495 ***
    // Wavefunction(s) for diagram number 13342
    // (none)
    // Amplitude(s) for diagram number 13342
    FFV1_0( w_fp[221], w_fp[694], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13343 OF 15495 ***
    // Wavefunction(s) for diagram number 13343
    // (none)
    // Amplitude(s) for diagram number 13343
    FFV1_0( w_fp[603], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13344 OF 15495 ***
    // Wavefunction(s) for diagram number 13344
    // (none)
    // Amplitude(s) for diagram number 13344
    VVV1_0( w_fp[136], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVV1_0( w_fp[137], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[551], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 13345 OF 15495 ***
    // Wavefunction(s) for diagram number 13345
    // (none)
    // Amplitude(s) for diagram number 13345
    FFV1_0( w_fp[3], w_fp[227], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13346 OF 15495 ***
    // Wavefunction(s) for diagram number 13346
    // (none)
    // Amplitude(s) for diagram number 13346
    FFV1_0( w_fp[3], w_fp[473], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13347 OF 15495 ***
    // Wavefunction(s) for diagram number 13347
    // (none)
    // Amplitude(s) for diagram number 13347
    VVV1_0( w_fp[519], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[454], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    VVV1_0( w_fp[575], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 13348 OF 15495 ***
    // Wavefunction(s) for diagram number 13348
    // (none)
    // Amplitude(s) for diagram number 13348
    VVV1_0( w_fp[0], w_fp[247], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[251], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[250], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 13349 OF 15495 ***
    // Wavefunction(s) for diagram number 13349
    // (none)
    // Amplitude(s) for diagram number 13349
    VVV1_0( w_fp[362], w_fp[5], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13350 OF 15495 ***
    // Wavefunction(s) for diagram number 13350
    // (none)
    // Amplitude(s) for diagram number 13350
    FFV1_0( w_fp[3], w_fp[690], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];

    // *** DIAGRAM 13351 OF 15495 ***
    // Wavefunction(s) for diagram number 13351
    // (none)
    // Amplitude(s) for diagram number 13351
    FFV1_0( w_fp[206], w_fp[75], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13352 OF 15495 ***
    // Wavefunction(s) for diagram number 13352
    // (none)
    // Amplitude(s) for diagram number 13352
    FFV1_0( w_fp[206], w_fp[690], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13353 OF 15495 ***
    // Wavefunction(s) for diagram number 13353
    // (none)
    // Amplitude(s) for diagram number 13353
    FFV1_0( w_fp[3], w_fp[75], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];

    // *** DIAGRAM 13354 OF 15495 ***
    // Wavefunction(s) for diagram number 13354
    // (none)
    // Amplitude(s) for diagram number 13354
    VVV1_0( w_fp[1], w_fp[87], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13355 OF 15495 ***
    // Wavefunction(s) for diagram number 13355
    // (none)
    // Amplitude(s) for diagram number 13355
    FFV1_0( w_fp[3], w_fp[542], w_fp[82], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[73], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[366], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13356 OF 15495 ***
    // Wavefunction(s) for diagram number 13356
    // (none)
    // Amplitude(s) for diagram number 13356
    FFV1_0( w_fp[607], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];

    // *** DIAGRAM 13357 OF 15495 ***
    // Wavefunction(s) for diagram number 13357
    // (none)
    // Amplitude(s) for diagram number 13357
    FFV1_0( w_fp[3], w_fp[473], w_fp[605], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13358 OF 15495 ***
    // Wavefunction(s) for diagram number 13358
    // (none)
    // Amplitude(s) for diagram number 13358
    VVVV1_0( w_fp[576], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVVV3_0( w_fp[576], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVVV4_0( w_fp[576], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 13359 OF 15495 ***
    // Wavefunction(s) for diagram number 13359
    // (none)
    // Amplitude(s) for diagram number 13359
    VVV1_0( w_fp[228], w_fp[5], w_fp[604], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13360 OF 15495 ***
    // Wavefunction(s) for diagram number 13360
    // (none)
    // Amplitude(s) for diagram number 13360
    VVV1_0( w_fp[1], w_fp[228], w_fp[605], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 13361 OF 15495 ***
    // Wavefunction(s) for diagram number 13361
    // (none)
    // Amplitude(s) for diagram number 13361
    FFV1_0( w_fp[3], w_fp[226], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13362 OF 15495 ***
    // Wavefunction(s) for diagram number 13362
    // (none)
    // Amplitude(s) for diagram number 13362
    FFV1_0( w_fp[607], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 13363 OF 15495 ***
    // Wavefunction(s) for diagram number 13363
    // (none)
    // Amplitude(s) for diagram number 13363
    FFV1_0( w_fp[206], w_fp[192], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13364 OF 15495 ***
    // Wavefunction(s) for diagram number 13364
    // (none)
    // Amplitude(s) for diagram number 13364
    FFV1_0( w_fp[610], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13365 OF 15495 ***
    // Wavefunction(s) for diagram number 13365
    // (none)
    // Amplitude(s) for diagram number 13365
    FFV1_0( w_fp[3], w_fp[192], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 13366 OF 15495 ***
    // Wavefunction(s) for diagram number 13366
    // (none)
    // Amplitude(s) for diagram number 13366
    FFV1_0( w_fp[3], w_fp[473], w_fp[608], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13367 OF 15495 ***
    // Wavefunction(s) for diagram number 13367
    // (none)
    // Amplitude(s) for diagram number 13367
    VVVV1_0( w_fp[0], w_fp[362], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[362], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[362], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 13368 OF 15495 ***
    // Wavefunction(s) for diagram number 13368
    // (none)
    // Amplitude(s) for diagram number 13368
    VVV1_0( w_fp[228], w_fp[5], w_fp[64], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 13369 OF 15495 ***
    // Wavefunction(s) for diagram number 13369
    // (none)
    // Amplitude(s) for diagram number 13369
    VVV1_0( w_fp[362], w_fp[5], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 13370 OF 15495 ***
    // Wavefunction(s) for diagram number 13370
    // (none)
    // Amplitude(s) for diagram number 13370
    FFV1_0( w_fp[3], w_fp[226], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13371 OF 15495 ***
    // Wavefunction(s) for diagram number 13371
    // (none)
    // Amplitude(s) for diagram number 13371
    FFV1_0( w_fp[3], w_fp[693], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 13372 OF 15495 ***
    // Wavefunction(s) for diagram number 13372
    // (none)
    // Amplitude(s) for diagram number 13372
    VVVV1_0( w_fp[0], w_fp[1], w_fp[228], w_fp[87], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[228], w_fp[87], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[228], w_fp[87], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 13373 OF 15495 ***
    // Wavefunction(s) for diagram number 13373
    // (none)
    // Amplitude(s) for diagram number 13373
    VVV1_0( w_fp[1], w_fp[87], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 13374 OF 15495 ***
    // Wavefunction(s) for diagram number 13374
    // (none)
    // Amplitude(s) for diagram number 13374
    VVV1_0( w_fp[1], w_fp[228], w_fp[608], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 13375 OF 15495 ***
    // Wavefunction(s) for diagram number 13375
    // (none)
    // Amplitude(s) for diagram number 13375
    FFV1_0( w_fp[206], w_fp[693], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13376 OF 15495 ***
    // Wavefunction(s) for diagram number 13376
    // (none)
    // Amplitude(s) for diagram number 13376
    FFV1_0( w_fp[610], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13377 OF 15495 ***
    // Wavefunction(s) for diagram number 13377
    // (none)
    // Amplitude(s) for diagram number 13377
    VVV1_0( w_fp[22], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0( w_fp[21], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVV1_0( w_fp[20], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13378 OF 15495 ***
    // Wavefunction(s) for diagram number 13378
    // (none)
    // Amplitude(s) for diagram number 13378
    FFV1_0( w_fp[3], w_fp[226], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13379 OF 15495 ***
    // Wavefunction(s) for diagram number 13379
    // (none)
    // Amplitude(s) for diagram number 13379
    FFV1_0( w_fp[3], w_fp[473], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[447], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[590], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13380 OF 15495 ***
    // Wavefunction(s) for diagram number 13380
    // (none)
    // Amplitude(s) for diagram number 13380
    VVV1_0( w_fp[561], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVV1_0( w_fp[447], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVV1_0( w_fp[590], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 13381 OF 15495 ***
    // Wavefunction(s) for diagram number 13381
    // (none)
    // Amplitude(s) for diagram number 13381
    VVV1_0( w_fp[0], w_fp[82], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[73], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[366], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 13382 OF 15495 ***
    // Wavefunction(s) for diagram number 13382
    // (none)
    // Amplitude(s) for diagram number 13382
    VVV1_0( w_fp[359], w_fp[4], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13383 OF 15495 ***
    // Wavefunction(s) for diagram number 13383
    // (none)
    // Amplitude(s) for diagram number 13383
    FFV1_0( w_fp[3], w_fp[692], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];

    // *** DIAGRAM 13384 OF 15495 ***
    // Wavefunction(s) for diagram number 13384
    // (none)
    // Amplitude(s) for diagram number 13384
    FFV1_0( w_fp[184], w_fp[75], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13385 OF 15495 ***
    // Wavefunction(s) for diagram number 13385
    // (none)
    // Amplitude(s) for diagram number 13385
    FFV1_0( w_fp[184], w_fp[692], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13386 OF 15495 ***
    // Wavefunction(s) for diagram number 13386
    // (none)
    // Amplitude(s) for diagram number 13386
    FFV1_0( w_fp[3], w_fp[75], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 13387 OF 15495 ***
    // Wavefunction(s) for diagram number 13387
    // (none)
    // Amplitude(s) for diagram number 13387
    VVV1_0( w_fp[1], w_fp[115], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13388 OF 15495 ***
    // Wavefunction(s) for diagram number 13388
    // (none)
    // Amplitude(s) for diagram number 13388
    FFV1_0( w_fp[3], w_fp[542], w_fp[262], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[353], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[72], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13389 OF 15495 ***
    // Wavefunction(s) for diagram number 13389
    // (none)
    // Amplitude(s) for diagram number 13389
    FFV1_0( w_fp[620], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];

    // *** DIAGRAM 13390 OF 15495 ***
    // Wavefunction(s) for diagram number 13390
    // (none)
    // Amplitude(s) for diagram number 13390
    FFV1_0( w_fp[3], w_fp[473], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13391 OF 15495 ***
    // Wavefunction(s) for diagram number 13391
    // (none)
    // Amplitude(s) for diagram number 13391
    VVVV1_0( w_fp[559], w_fp[1], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVVV3_0( w_fp[559], w_fp[1], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVVV4_0( w_fp[559], w_fp[1], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 13392 OF 15495 ***
    // Wavefunction(s) for diagram number 13392
    // (none)
    // Amplitude(s) for diagram number 13392
    VVV1_0( w_fp[228], w_fp[4], w_fp[120], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13393 OF 15495 ***
    // Wavefunction(s) for diagram number 13393
    // (none)
    // Amplitude(s) for diagram number 13393
    VVV1_0( w_fp[1], w_fp[228], w_fp[585], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 13394 OF 15495 ***
    // Wavefunction(s) for diagram number 13394
    // (none)
    // Amplitude(s) for diagram number 13394
    FFV1_0( w_fp[3], w_fp[225], w_fp[120], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13395 OF 15495 ***
    // Wavefunction(s) for diagram number 13395
    // (none)
    // Amplitude(s) for diagram number 13395
    FFV1_0( w_fp[620], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 13396 OF 15495 ***
    // Wavefunction(s) for diagram number 13396
    // (none)
    // Amplitude(s) for diagram number 13396
    FFV1_0( w_fp[184], w_fp[192], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13397 OF 15495 ***
    // Wavefunction(s) for diagram number 13397
    // (none)
    // Amplitude(s) for diagram number 13397
    FFV1_0( w_fp[623], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13398 OF 15495 ***
    // Wavefunction(s) for diagram number 13398
    // (none)
    // Amplitude(s) for diagram number 13398
    FFV1_0( w_fp[3], w_fp[192], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 13399 OF 15495 ***
    // Wavefunction(s) for diagram number 13399
    // (none)
    // Amplitude(s) for diagram number 13399
    FFV1_0( w_fp[3], w_fp[473], w_fp[621], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13400 OF 15495 ***
    // Wavefunction(s) for diagram number 13400
    // (none)
    // Amplitude(s) for diagram number 13400
    VVVV1_0( w_fp[0], w_fp[359], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[359], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[359], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    storeWf( wfs, w_cx, nevt, 75 );
    storeWf( wfs, w_cx, nevt, 192 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
