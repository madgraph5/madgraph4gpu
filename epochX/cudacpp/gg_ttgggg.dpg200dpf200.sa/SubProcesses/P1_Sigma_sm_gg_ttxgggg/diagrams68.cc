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
  diagramgroup68( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 22 );
    retrieveWf( wfs, w_cx, nevt, 23 );
    retrieveWf( wfs, w_cx, nevt, 24 );
    retrieveWf( wfs, w_cx, nevt, 25 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 48 );
    retrieveWf( wfs, w_cx, nevt, 49 );
    retrieveWf( wfs, w_cx, nevt, 50 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 106 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 117 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 121 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 126 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 165 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 229 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 240 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 249 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 257 );
    retrieveWf( wfs, w_cx, nevt, 262 );
    retrieveWf( wfs, w_cx, nevt, 264 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 291 );
    retrieveWf( wfs, w_cx, nevt, 345 );
    retrieveWf( wfs, w_cx, nevt, 353 );
    retrieveWf( wfs, w_cx, nevt, 354 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 358 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 367 );
    retrieveWf( wfs, w_cx, nevt, 368 );
    retrieveWf( wfs, w_cx, nevt, 369 );
    retrieveWf( wfs, w_cx, nevt, 371 );
    retrieveWf( wfs, w_cx, nevt, 372 );
    retrieveWf( wfs, w_cx, nevt, 373 );
    retrieveWf( wfs, w_cx, nevt, 374 );
    retrieveWf( wfs, w_cx, nevt, 375 );
    retrieveWf( wfs, w_cx, nevt, 432 );
    retrieveWf( wfs, w_cx, nevt, 433 );
    retrieveWf( wfs, w_cx, nevt, 434 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 499 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 517 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 563 );
    retrieveWf( wfs, w_cx, nevt, 564 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 572 );
    retrieveWf( wfs, w_cx, nevt, 573 );
    retrieveWf( wfs, w_cx, nevt, 574 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 587 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 589 );
    retrieveWf( wfs, w_cx, nevt, 600 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 606 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 611 );
    retrieveWf( wfs, w_cx, nevt, 613 );
    retrieveWf( wfs, w_cx, nevt, 616 );
    retrieveWf( wfs, w_cx, nevt, 619 );
    retrieveWf( wfs, w_cx, nevt, 621 );
    retrieveWf( wfs, w_cx, nevt, 622 );
    retrieveWf( wfs, w_cx, nevt, 623 );
    retrieveWf( wfs, w_cx, nevt, 626 );
    retrieveWf( wfs, w_cx, nevt, 651 );
    retrieveWf( wfs, w_cx, nevt, 652 );
    retrieveWf( wfs, w_cx, nevt, 653 );
    retrieveWf( wfs, w_cx, nevt, 658 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 663 );
    retrieveWf( wfs, w_cx, nevt, 664 );
    retrieveWf( wfs, w_cx, nevt, 673 );
    retrieveWf( wfs, w_cx, nevt, 674 );
    retrieveWf( wfs, w_cx, nevt, 675 );
    retrieveWf( wfs, w_cx, nevt, 695 );
    retrieveWf( wfs, w_cx, nevt, 697 );
    retrieveWf( wfs, w_cx, nevt, 700 );
    retrieveWf( wfs, w_cx, nevt, 746 );
    retrieveWf( wfs, w_cx, nevt, 747 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 13401 OF 15495 ***
    // Wavefunction(s) for diagram number 13401
    // (none)
    // Amplitude(s) for diagram number 13401
    VVV1_0( w_fp[228], w_fp[4], w_fp[114], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13402 OF 15495 ***
    // Wavefunction(s) for diagram number 13402
    // (none)
    // Amplitude(s) for diagram number 13402
    VVV1_0( w_fp[359], w_fp[4], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13403 OF 15495 ***
    // Wavefunction(s) for diagram number 13403
    // (none)
    // Amplitude(s) for diagram number 13403
    FFV1_0( w_fp[3], w_fp[225], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13404 OF 15495 ***
    // Wavefunction(s) for diagram number 13404
    // (none)
    // Amplitude(s) for diagram number 13404
    FFV1_0( w_fp[3], w_fp[695], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 13405 OF 15495 ***
    // Wavefunction(s) for diagram number 13405
    // (none)
    // Amplitude(s) for diagram number 13405
    VVVV1_0( w_fp[0], w_fp[1], w_fp[228], w_fp[115], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[228], w_fp[115], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[228], w_fp[115], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 13406 OF 15495 ***
    // Wavefunction(s) for diagram number 13406
    // (none)
    // Amplitude(s) for diagram number 13406
    VVV1_0( w_fp[1], w_fp[115], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 13407 OF 15495 ***
    // Wavefunction(s) for diagram number 13407
    // (none)
    // Amplitude(s) for diagram number 13407
    VVV1_0( w_fp[1], w_fp[228], w_fp[621], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 13408 OF 15495 ***
    // Wavefunction(s) for diagram number 13408
    // (none)
    // Amplitude(s) for diagram number 13408
    FFV1_0( w_fp[184], w_fp[695], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13409 OF 15495 ***
    // Wavefunction(s) for diagram number 13409
    // (none)
    // Amplitude(s) for diagram number 13409
    FFV1_0( w_fp[623], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13410 OF 15495 ***
    // Wavefunction(s) for diagram number 13410
    // (none)
    // Amplitude(s) for diagram number 13410
    VVV1_0( w_fp[106], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[117], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    VVV1_0( w_fp[59], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13411 OF 15495 ***
    // Wavefunction(s) for diagram number 13411
    // (none)
    // Amplitude(s) for diagram number 13411
    FFV1_0( w_fp[3], w_fp[225], w_fp[106], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[225], w_fp[117], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[225], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13412 OF 15495 ***
    // Wavefunction(s) for diagram number 13412
    // (none)
    // Amplitude(s) for diagram number 13412
    FFV1_0( w_fp[3], w_fp[473], w_fp[517], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[473], w_fp[472], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[445], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13413 OF 15495 ***
    // Wavefunction(s) for diagram number 13413
    // (none)
    // Amplitude(s) for diagram number 13413
    VVV1_0( w_fp[517], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[472], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    VVV1_0( w_fp[445], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13414 OF 15495 ***
    // Wavefunction(s) for diagram number 13414
    // (none)
    // Amplitude(s) for diagram number 13414
    VVV1_0( w_fp[0], w_fp[262], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[353], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[0], w_fp[72], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 13415 OF 15495 ***
    // Wavefunction(s) for diagram number 13415
    // (none)
    // Amplitude(s) for diagram number 13415
    FFV1_0( w_fp[3], w_fp[542], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[542], w_fp[433], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[432], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13416 OF 15495 ***
    // Wavefunction(s) for diagram number 13416
    // (none)
    // Amplitude(s) for diagram number 13416
    FFV1_0( w_fp[50], w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    FFV1_0( w_fp[49], w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    FFV1_0( w_fp[48], w_fp[542], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];

    // *** DIAGRAM 13417 OF 15495 ***
    // Wavefunction(s) for diagram number 13417
    // (none)
    // Amplitude(s) for diagram number 13417
    FFV1_0( w_fp[3], w_fp[473], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[473], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13418 OF 15495 ***
    // Wavefunction(s) for diagram number 13418
    // (none)
    // Amplitude(s) for diagram number 13418
    VVV1_0( w_fp[515], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[589], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    VVV1_0( w_fp[510], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 13419 OF 15495 ***
    // Wavefunction(s) for diagram number 13419
    // (none)
    // Amplitude(s) for diagram number 13419
    FFV1_0( w_fp[50], w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    FFV1_0( w_fp[49], w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    FFV1_0( w_fp[48], w_fp[473], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 13420 OF 15495 ***
    // Wavefunction(s) for diagram number 13420
    // (none)
    // Amplitude(s) for diagram number 13420
    VVV1_0( w_fp[0], w_fp[91], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[433], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[0], w_fp[432], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13421 OF 15495 ***
    // Wavefunction(s) for diagram number 13421
    // (none)
    // Amplitude(s) for diagram number 13421
    FFV1_0( w_fp[3], w_fp[215], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[240], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[165], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[574], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[573], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[572], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[563], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[564], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 13422 OF 15495 ***
    // Wavefunction(s) for diagram number 13422
    // (none)
    // Amplitude(s) for diagram number 13422
    VVV1_0( w_fp[359], w_fp[7], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13423 OF 15495 ***
    // Wavefunction(s) for diagram number 13423
    // (none)
    // Amplitude(s) for diagram number 13423
    FFV1_0( w_fp[45], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];

    // *** DIAGRAM 13424 OF 15495 ***
    // Wavefunction(s) for diagram number 13424
    // (none)
    // Amplitude(s) for diagram number 13424
    FFV1_0( w_fp[160], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13425 OF 15495 ***
    // Wavefunction(s) for diagram number 13425
    // (none)
    // Amplitude(s) for diagram number 13425
    FFV1_0( w_fp[45], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13426 OF 15495 ***
    // Wavefunction(s) for diagram number 13426
    // (none)
    // Amplitude(s) for diagram number 13426
    FFV1_0( w_fp[160], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13427 OF 15495 ***
    // Wavefunction(s) for diagram number 13427
    // (none)
    // Amplitude(s) for diagram number 13427
    VVV1_0( w_fp[1], w_fp[116], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13428 OF 15495 ***
    // Wavefunction(s) for diagram number 13428
    // (none)
    // Amplitude(s) for diagram number 13428
    FFV1_0( w_fp[536], w_fp[2], w_fp[79], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[374], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[375], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13429 OF 15495 ***
    // Wavefunction(s) for diagram number 13429
    // (none)
    // Amplitude(s) for diagram number 13429
    FFV1_0( w_fp[355], w_fp[664], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];

    // *** DIAGRAM 13430 OF 15495 ***
    // Wavefunction(s) for diagram number 13430
    // (none)
    // Amplitude(s) for diagram number 13430
    FFV1_0( w_fp[355], w_fp[2], w_fp[619], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13431 OF 15495 ***
    // Wavefunction(s) for diagram number 13431
    // (none)
    // Amplitude(s) for diagram number 13431
    VVVV1_0( w_fp[559], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    VVVV3_0( w_fp[559], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    VVVV4_0( w_fp[559], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13432 OF 15495 ***
    // Wavefunction(s) for diagram number 13432
    // (none)
    // Amplitude(s) for diagram number 13432
    VVV1_0( w_fp[239], w_fp[7], w_fp[120], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13433 OF 15495 ***
    // Wavefunction(s) for diagram number 13433
    // (none)
    // Amplitude(s) for diagram number 13433
    VVV1_0( w_fp[1], w_fp[239], w_fp[619], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13434 OF 15495 ***
    // Wavefunction(s) for diagram number 13434
    // (none)
    // Amplitude(s) for diagram number 13434
    FFV1_0( w_fp[199], w_fp[2], w_fp[120], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13435 OF 15495 ***
    // Wavefunction(s) for diagram number 13435
    // (none)
    // Amplitude(s) for diagram number 13435
    FFV1_0( w_fp[199], w_fp[664], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];

    // *** DIAGRAM 13436 OF 15495 ***
    // Wavefunction(s) for diagram number 13436
    // (none)
    // Amplitude(s) for diagram number 13436
    FFV1_0( w_fp[674], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13437 OF 15495 ***
    // Wavefunction(s) for diagram number 13437
    // (none)
    // Amplitude(s) for diagram number 13437
    FFV1_0( w_fp[355], w_fp[121], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13438 OF 15495 ***
    // Wavefunction(s) for diagram number 13438
    // (none)
    // Amplitude(s) for diagram number 13438
    FFV1_0( w_fp[674], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13439 OF 15495 ***
    // Wavefunction(s) for diagram number 13439
    // (none)
    // Amplitude(s) for diagram number 13439
    FFV1_0( w_fp[355], w_fp[2], w_fp[622], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13440 OF 15495 ***
    // Wavefunction(s) for diagram number 13440
    // (none)
    // Amplitude(s) for diagram number 13440
    VVVV1_0( w_fp[0], w_fp[359], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[359], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[359], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 13441 OF 15495 ***
    // Wavefunction(s) for diagram number 13441
    // (none)
    // Amplitude(s) for diagram number 13441
    VVV1_0( w_fp[239], w_fp[7], w_fp[114], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13442 OF 15495 ***
    // Wavefunction(s) for diagram number 13442
    // (none)
    // Amplitude(s) for diagram number 13442
    VVV1_0( w_fp[359], w_fp[7], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13443 OF 15495 ***
    // Wavefunction(s) for diagram number 13443
    // (none)
    // Amplitude(s) for diagram number 13443
    FFV1_0( w_fp[199], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13444 OF 15495 ***
    // Wavefunction(s) for diagram number 13444
    // (none)
    // Amplitude(s) for diagram number 13444
    FFV1_0( w_fp[532], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];

    // *** DIAGRAM 13445 OF 15495 ***
    // Wavefunction(s) for diagram number 13445
    // (none)
    // Amplitude(s) for diagram number 13445
    VVVV1_0( w_fp[0], w_fp[1], w_fp[239], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[239], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[239], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13446 OF 15495 ***
    // Wavefunction(s) for diagram number 13446
    // (none)
    // Amplitude(s) for diagram number 13446
    VVV1_0( w_fp[1], w_fp[116], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13447 OF 15495 ***
    // Wavefunction(s) for diagram number 13447
    // (none)
    // Amplitude(s) for diagram number 13447
    VVV1_0( w_fp[1], w_fp[239], w_fp[622], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13448 OF 15495 ***
    // Wavefunction(s) for diagram number 13448
    // (none)
    // Amplitude(s) for diagram number 13448
    FFV1_0( w_fp[199], w_fp[121], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13449 OF 15495 ***
    // Wavefunction(s) for diagram number 13449
    // (none)
    // Amplitude(s) for diagram number 13449
    FFV1_0( w_fp[532], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13450 OF 15495 ***
    // Wavefunction(s) for diagram number 13450
    // (none)
    // Amplitude(s) for diagram number 13450
    VVV1_0( w_fp[106], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVV1_0( w_fp[117], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0( w_fp[59], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 13451 OF 15495 ***
    // Wavefunction(s) for diagram number 13451
    // (none)
    // Amplitude(s) for diagram number 13451
    FFV1_0( w_fp[199], w_fp[2], w_fp[106], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[117], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13452 OF 15495 ***
    // Wavefunction(s) for diagram number 13452
    // (none)
    // Amplitude(s) for diagram number 13452
    FFV1_0( w_fp[355], w_fp[2], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13453 OF 15495 ***
    // Wavefunction(s) for diagram number 13453
    // (none)
    // Amplitude(s) for diagram number 13453
    VVV1_0( w_fp[479], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0( w_fp[474], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    VVV1_0( w_fp[444], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13454 OF 15495 ***
    // Wavefunction(s) for diagram number 13454
    // (none)
    // Amplitude(s) for diagram number 13454
    VVV1_0( w_fp[0], w_fp[79], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[374], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[375], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13455 OF 15495 ***
    // Wavefunction(s) for diagram number 13455
    // (none)
    // Amplitude(s) for diagram number 13455
    VVV1_0( w_fp[360], w_fp[6], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13456 OF 15495 ***
    // Wavefunction(s) for diagram number 13456
    // (none)
    // Amplitude(s) for diagram number 13456
    FFV1_0( w_fp[560], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];

    // *** DIAGRAM 13457 OF 15495 ***
    // Wavefunction(s) for diagram number 13457
    // (none)
    // Amplitude(s) for diagram number 13457
    FFV1_0( w_fp[160], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13458 OF 15495 ***
    // Wavefunction(s) for diagram number 13458
    // (none)
    // Amplitude(s) for diagram number 13458
    FFV1_0( w_fp[560], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13459 OF 15495 ***
    // Wavefunction(s) for diagram number 13459
    // (none)
    // Amplitude(s) for diagram number 13459
    FFV1_0( w_fp[160], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];

    // *** DIAGRAM 13460 OF 15495 ***
    // Wavefunction(s) for diagram number 13460
    // (none)
    // Amplitude(s) for diagram number 13460
    VVV1_0( w_fp[1], w_fp[125], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13461 OF 15495 ***
    // Wavefunction(s) for diagram number 13461
    // (none)
    // Amplitude(s) for diagram number 13461
    FFV1_0( w_fp[536], w_fp[2], w_fp[358], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[81], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13462 OF 15495 ***
    // Wavefunction(s) for diagram number 13462
    // (none)
    // Amplitude(s) for diagram number 13462
    FFV1_0( w_fp[355], w_fp[663], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];

    // *** DIAGRAM 13463 OF 15495 ***
    // Wavefunction(s) for diagram number 13463
    // (none)
    // Amplitude(s) for diagram number 13463
    FFV1_0( w_fp[355], w_fp[2], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13464 OF 15495 ***
    // Wavefunction(s) for diagram number 13464
    // (none)
    // Amplitude(s) for diagram number 13464
    VVVV1_0( w_fp[522], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    VVVV3_0( w_fp[522], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    VVVV4_0( w_fp[522], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];

    // *** DIAGRAM 13465 OF 15495 ***
    // Wavefunction(s) for diagram number 13465
    // (none)
    // Amplitude(s) for diagram number 13465
    VVV1_0( w_fp[239], w_fp[6], w_fp[673], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 13466 OF 15495 ***
    // Wavefunction(s) for diagram number 13466
    // (none)
    // Amplitude(s) for diagram number 13466
    VVV1_0( w_fp[1], w_fp[239], w_fp[505], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];

    // *** DIAGRAM 13467 OF 15495 ***
    // Wavefunction(s) for diagram number 13467
    // (none)
    // Amplitude(s) for diagram number 13467
    FFV1_0( w_fp[198], w_fp[2], w_fp[673], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13468 OF 15495 ***
    // Wavefunction(s) for diagram number 13468
    // (none)
    // Amplitude(s) for diagram number 13468
    FFV1_0( w_fp[198], w_fp[663], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 13469 OF 15495 ***
    // Wavefunction(s) for diagram number 13469
    // (none)
    // Amplitude(s) for diagram number 13469
    FFV1_0( w_fp[674], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13470 OF 15495 ***
    // Wavefunction(s) for diagram number 13470
    // (none)
    // Amplitude(s) for diagram number 13470
    FFV1_0( w_fp[355], w_fp[675], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13471 OF 15495 ***
    // Wavefunction(s) for diagram number 13471
    // (none)
    // Amplitude(s) for diagram number 13471
    FFV1_0( w_fp[674], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];

    // *** DIAGRAM 13472 OF 15495 ***
    // Wavefunction(s) for diagram number 13472
    // (none)
    // Amplitude(s) for diagram number 13472
    FFV1_0( w_fp[355], w_fp[2], w_fp[626], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13473 OF 15495 ***
    // Wavefunction(s) for diagram number 13473
    // (none)
    // Amplitude(s) for diagram number 13473
    VVVV1_0( w_fp[0], w_fp[360], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[360], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[360], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];

    // *** DIAGRAM 13474 OF 15495 ***
    // Wavefunction(s) for diagram number 13474
    // (none)
    // Amplitude(s) for diagram number 13474
    VVV1_0( w_fp[239], w_fp[6], w_fp[55], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 13475 OF 15495 ***
    // Wavefunction(s) for diagram number 13475
    // (none)
    // Amplitude(s) for diagram number 13475
    VVV1_0( w_fp[360], w_fp[6], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];

    // *** DIAGRAM 13476 OF 15495 ***
    // Wavefunction(s) for diagram number 13476
    // (none)
    // Amplitude(s) for diagram number 13476
    FFV1_0( w_fp[198], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13477 OF 15495 ***
    // Wavefunction(s) for diagram number 13477
    // (none)
    // Amplitude(s) for diagram number 13477
    FFV1_0( w_fp[541], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];

    // *** DIAGRAM 13478 OF 15495 ***
    // Wavefunction(s) for diagram number 13478
    // (none)
    // Amplitude(s) for diagram number 13478
    VVVV1_0( w_fp[0], w_fp[1], w_fp[239], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[239], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[239], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];

    // *** DIAGRAM 13479 OF 15495 ***
    // Wavefunction(s) for diagram number 13479
    // (none)
    // Amplitude(s) for diagram number 13479
    VVV1_0( w_fp[1], w_fp[125], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];

    // *** DIAGRAM 13480 OF 15495 ***
    // Wavefunction(s) for diagram number 13480
    // (none)
    // Amplitude(s) for diagram number 13480
    VVV1_0( w_fp[1], w_fp[239], w_fp[626], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];

    // *** DIAGRAM 13481 OF 15495 ***
    // Wavefunction(s) for diagram number 13481
    // (none)
    // Amplitude(s) for diagram number 13481
    FFV1_0( w_fp[198], w_fp[675], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13482 OF 15495 ***
    // Wavefunction(s) for diagram number 13482
    // (none)
    // Amplitude(s) for diagram number 13482
    FFV1_0( w_fp[541], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13483 OF 15495 ***
    // Wavefunction(s) for diagram number 13483
    // (none)
    // Amplitude(s) for diagram number 13483
    VVV1_0( w_fp[90], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    VVV1_0( w_fp[126], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    VVV1_0( w_fp[62], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];

    // *** DIAGRAM 13484 OF 15495 ***
    // Wavefunction(s) for diagram number 13484
    // (none)
    // Amplitude(s) for diagram number 13484
    FFV1_0( w_fp[198], w_fp[2], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13485 OF 15495 ***
    // Wavefunction(s) for diagram number 13485
    // (none)
    // Amplitude(s) for diagram number 13485
    FFV1_0( w_fp[355], w_fp[2], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13486 OF 15495 ***
    // Wavefunction(s) for diagram number 13486
    // (none)
    // Amplitude(s) for diagram number 13486
    VVV1_0( w_fp[523], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    VVV1_0( w_fp[580], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    VVV1_0( w_fp[528], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];

    // *** DIAGRAM 13487 OF 15495 ***
    // Wavefunction(s) for diagram number 13487
    // (none)
    // Amplitude(s) for diagram number 13487
    VVV1_0( w_fp[0], w_fp[358], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[81], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[46], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];

    // *** DIAGRAM 13488 OF 15495 ***
    // Wavefunction(s) for diagram number 13488
    // (none)
    // Amplitude(s) for diagram number 13488
    VVV1_0( w_fp[361], w_fp[5], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13489 OF 15495 ***
    // Wavefunction(s) for diagram number 13489
    // (none)
    // Amplitude(s) for diagram number 13489
    FFV1_0( w_fp[586], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];

    // *** DIAGRAM 13490 OF 15495 ***
    // Wavefunction(s) for diagram number 13490
    // (none)
    // Amplitude(s) for diagram number 13490
    FFV1_0( w_fp[160], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13491 OF 15495 ***
    // Wavefunction(s) for diagram number 13491
    // (none)
    // Amplitude(s) for diagram number 13491
    FFV1_0( w_fp[586], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13492 OF 15495 ***
    // Wavefunction(s) for diagram number 13492
    // (none)
    // Amplitude(s) for diagram number 13492
    FFV1_0( w_fp[160], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13493 OF 15495 ***
    // Wavefunction(s) for diagram number 13493
    // (none)
    // Amplitude(s) for diagram number 13493
    VVV1_0( w_fp[1], w_fp[131], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13494 OF 15495 ***
    // Wavefunction(s) for diagram number 13494
    // (none)
    // Amplitude(s) for diagram number 13494
    FFV1_0( w_fp[536], w_fp[2], w_fp[257], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[249], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[354], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13495 OF 15495 ***
    // Wavefunction(s) for diagram number 13495
    // (none)
    // Amplitude(s) for diagram number 13495
    FFV1_0( w_fp[355], w_fp[497], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 13496 OF 15495 ***
    // Wavefunction(s) for diagram number 13496
    // (none)
    // Amplitude(s) for diagram number 13496
    FFV1_0( w_fp[355], w_fp[2], w_fp[587], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13497 OF 15495 ***
    // Wavefunction(s) for diagram number 13497
    // (none)
    // Amplitude(s) for diagram number 13497
    VVVV1_0( w_fp[516], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    VVVV3_0( w_fp[516], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    VVVV4_0( w_fp[516], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 13498 OF 15495 ***
    // Wavefunction(s) for diagram number 13498
    // (none)
    // Amplitude(s) for diagram number 13498
    VVV1_0( w_fp[239], w_fp[5], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 13499 OF 15495 ***
    // Wavefunction(s) for diagram number 13499
    // (none)
    // Amplitude(s) for diagram number 13499
    VVV1_0( w_fp[1], w_fp[239], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 13500 OF 15495 ***
    // Wavefunction(s) for diagram number 13500
    // (none)
    // Amplitude(s) for diagram number 13500
    FFV1_0( w_fp[216], w_fp[2], w_fp[700], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13501 OF 15495 ***
    // Wavefunction(s) for diagram number 13501
    // (none)
    // Amplitude(s) for diagram number 13501
    FFV1_0( w_fp[216], w_fp[497], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 13502 OF 15495 ***
    // Wavefunction(s) for diagram number 13502
    // (none)
    // Amplitude(s) for diagram number 13502
    FFV1_0( w_fp[674], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13503 OF 15495 ***
    // Wavefunction(s) for diagram number 13503
    // (none)
    // Amplitude(s) for diagram number 13503
    FFV1_0( w_fp[355], w_fp[697], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13504 OF 15495 ***
    // Wavefunction(s) for diagram number 13504
    // (none)
    // Amplitude(s) for diagram number 13504
    FFV1_0( w_fp[674], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13505 OF 15495 ***
    // Wavefunction(s) for diagram number 13505
    // (none)
    // Amplitude(s) for diagram number 13505
    FFV1_0( w_fp[355], w_fp[2], w_fp[600], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13506 OF 15495 ***
    // Wavefunction(s) for diagram number 13506
    // (none)
    // Amplitude(s) for diagram number 13506
    VVVV1_0( w_fp[0], w_fp[361], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[361], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[361], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 13507 OF 15495 ***
    // Wavefunction(s) for diagram number 13507
    // (none)
    // Amplitude(s) for diagram number 13507
    VVV1_0( w_fp[239], w_fp[5], w_fp[238], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 13508 OF 15495 ***
    // Wavefunction(s) for diagram number 13508
    // (none)
    // Amplitude(s) for diagram number 13508
    VVV1_0( w_fp[361], w_fp[5], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 13509 OF 15495 ***
    // Wavefunction(s) for diagram number 13509
    // (none)
    // Amplitude(s) for diagram number 13509
    FFV1_0( w_fp[216], w_fp[2], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13510 OF 15495 ***
    // Wavefunction(s) for diagram number 13510
    // (none)
    // Amplitude(s) for diagram number 13510
    FFV1_0( w_fp[588], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 13511 OF 15495 ***
    // Wavefunction(s) for diagram number 13511
    // (none)
    // Amplitude(s) for diagram number 13511
    VVVV1_0( w_fp[0], w_fp[1], w_fp[239], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[239], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[239], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13512 OF 15495 ***
    // Wavefunction(s) for diagram number 13512
    // (none)
    // Amplitude(s) for diagram number 13512
    VVV1_0( w_fp[1], w_fp[131], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13513 OF 15495 ***
    // Wavefunction(s) for diagram number 13513
    // (none)
    // Amplitude(s) for diagram number 13513
    VVV1_0( w_fp[1], w_fp[239], w_fp[600], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 13514 OF 15495 ***
    // Wavefunction(s) for diagram number 13514
    // (none)
    // Amplitude(s) for diagram number 13514
    FFV1_0( w_fp[216], w_fp[697], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13515 OF 15495 ***
    // Wavefunction(s) for diagram number 13515
    // (none)
    // Amplitude(s) for diagram number 13515
    FFV1_0( w_fp[588], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13516 OF 15495 ***
    // Wavefunction(s) for diagram number 13516
    // (none)
    // Amplitude(s) for diagram number 13516
    VVV1_0( w_fp[71], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVV1_0( w_fp[132], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    VVV1_0( w_fp[65], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 13517 OF 15495 ***
    // Wavefunction(s) for diagram number 13517
    // (none)
    // Amplitude(s) for diagram number 13517
    FFV1_0( w_fp[216], w_fp[2], w_fp[71], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[132], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13518 OF 15495 ***
    // Wavefunction(s) for diagram number 13518
    // (none)
    // Amplitude(s) for diagram number 13518
    FFV1_0( w_fp[355], w_fp[2], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13519 OF 15495 ***
    // Wavefunction(s) for diagram number 13519
    // (none)
    // Amplitude(s) for diagram number 13519
    VVV1_0( w_fp[549], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0( w_fp[451], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    VVV1_0( w_fp[571], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13520 OF 15495 ***
    // Wavefunction(s) for diagram number 13520
    // (none)
    // Amplitude(s) for diagram number 13520
    VVV1_0( w_fp[0], w_fp[257], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[249], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[354], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13521 OF 15495 ***
    // Wavefunction(s) for diagram number 13521
    // (none)
    // Amplitude(s) for diagram number 13521
    FFV1_0( w_fp[536], w_fp[2], w_fp[345], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[95], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[434], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13522 OF 15495 ***
    // Wavefunction(s) for diagram number 13522
    // (none)
    // Amplitude(s) for diagram number 13522
    FFV1_0( w_fp[536], w_fp[229], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    FFV1_0( w_fp[536], w_fp[28], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    FFV1_0( w_fp[536], w_fp[60], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 13523 OF 15495 ***
    // Wavefunction(s) for diagram number 13523
    // (none)
    // Amplitude(s) for diagram number 13523
    FFV1_0( w_fp[355], w_fp[2], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13524 OF 15495 ***
    // Wavefunction(s) for diagram number 13524
    // (none)
    // Amplitude(s) for diagram number 13524
    VVV1_0( w_fp[452], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0( w_fp[488], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    VVV1_0( w_fp[437], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13525 OF 15495 ***
    // Wavefunction(s) for diagram number 13525
    // (none)
    // Amplitude(s) for diagram number 13525
    FFV1_0( w_fp[355], w_fp[229], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    FFV1_0( w_fp[355], w_fp[28], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    FFV1_0( w_fp[355], w_fp[60], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 13526 OF 15495 ***
    // Wavefunction(s) for diagram number 13526
    // (none)
    // Amplitude(s) for diagram number 13526
    VVV1_0( w_fp[0], w_fp[345], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[95], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[434], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 13527 OF 15495 ***
    // Wavefunction(s) for diagram number 13527
    // (none)
    // Amplitude(s) for diagram number 13527
    FFV1_0( w_fp[196], w_fp[2], w_fp[23], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[24], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[25], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[653], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[652], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[651], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 13528 OF 15495 ***
    // Wavefunction(s) for diagram number 13528
    // (none)
    // Amplitude(s) for diagram number 13528
    VVV1_0( w_fp[362], w_fp[7], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13529 OF 15495 ***
    // Wavefunction(s) for diagram number 13529
    // (none)
    // Amplitude(s) for diagram number 13529
    FFV1_0( w_fp[543], w_fp[2], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];

    // *** DIAGRAM 13530 OF 15495 ***
    // Wavefunction(s) for diagram number 13530
    // (none)
    // Amplitude(s) for diagram number 13530
    FFV1_0( w_fp[660], w_fp[118], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13531 OF 15495 ***
    // Wavefunction(s) for diagram number 13531
    // (none)
    // Amplitude(s) for diagram number 13531
    FFV1_0( w_fp[543], w_fp[118], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13532 OF 15495 ***
    // Wavefunction(s) for diagram number 13532
    // (none)
    // Amplitude(s) for diagram number 13532
    FFV1_0( w_fp[660], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13533 OF 15495 ***
    // Wavefunction(s) for diagram number 13533
    // (none)
    // Amplitude(s) for diagram number 13533
    VVV1_0( w_fp[1], w_fp[88], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13534 OF 15495 ***
    // Wavefunction(s) for diagram number 13534
    // (none)
    // Amplitude(s) for diagram number 13534
    FFV1_0( w_fp[527], w_fp[2], w_fp[367], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[368], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[369], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13535 OF 15495 ***
    // Wavefunction(s) for diagram number 13535
    // (none)
    // Amplitude(s) for diagram number 13535
    FFV1_0( w_fp[256], w_fp[246], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];

    // *** DIAGRAM 13536 OF 15495 ***
    // Wavefunction(s) for diagram number 13536
    // (none)
    // Amplitude(s) for diagram number 13536
    FFV1_0( w_fp[256], w_fp[2], w_fp[606], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13537 OF 15495 ***
    // Wavefunction(s) for diagram number 13537
    // (none)
    // Amplitude(s) for diagram number 13537
    VVVV1_0( w_fp[576], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    VVVV3_0( w_fp[576], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    VVVV4_0( w_fp[576], w_fp[1], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13538 OF 15495 ***
    // Wavefunction(s) for diagram number 13538
    // (none)
    // Amplitude(s) for diagram number 13538
    VVV1_0( w_fp[154], w_fp[7], w_fp[604], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13539 OF 15495 ***
    // Wavefunction(s) for diagram number 13539
    // (none)
    // Amplitude(s) for diagram number 13539
    VVV1_0( w_fp[1], w_fp[154], w_fp[606], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13540 OF 15495 ***
    // Wavefunction(s) for diagram number 13540
    // (none)
    // Amplitude(s) for diagram number 13540
    FFV1_0( w_fp[171], w_fp[2], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13541 OF 15495 ***
    // Wavefunction(s) for diagram number 13541
    // (none)
    // Amplitude(s) for diagram number 13541
    FFV1_0( w_fp[171], w_fp[246], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 13542 OF 15495 ***
    // Wavefunction(s) for diagram number 13542
    // (none)
    // Amplitude(s) for diagram number 13542
    FFV1_0( w_fp[658], w_fp[118], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13543 OF 15495 ***
    // Wavefunction(s) for diagram number 13543
    // (none)
    // Amplitude(s) for diagram number 13543
    FFV1_0( w_fp[256], w_fp[92], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13544 OF 15495 ***
    // Wavefunction(s) for diagram number 13544
    // (none)
    // Amplitude(s) for diagram number 13544
    FFV1_0( w_fp[658], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 13545 OF 15495 ***
    // Wavefunction(s) for diagram number 13545
    // (none)
    // Amplitude(s) for diagram number 13545
    FFV1_0( w_fp[256], w_fp[2], w_fp[609], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13546 OF 15495 ***
    // Wavefunction(s) for diagram number 13546
    // (none)
    // Amplitude(s) for diagram number 13546
    VVVV1_0( w_fp[0], w_fp[362], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[362], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[362], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 13547 OF 15495 ***
    // Wavefunction(s) for diagram number 13547
    // (none)
    // Amplitude(s) for diagram number 13547
    VVV1_0( w_fp[154], w_fp[7], w_fp[64], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13548 OF 15495 ***
    // Wavefunction(s) for diagram number 13548
    // (none)
    // Amplitude(s) for diagram number 13548
    VVV1_0( w_fp[362], w_fp[7], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13549 OF 15495 ***
    // Wavefunction(s) for diagram number 13549
    // (none)
    // Amplitude(s) for diagram number 13549
    FFV1_0( w_fp[171], w_fp[2], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13550 OF 15495 ***
    // Wavefunction(s) for diagram number 13550
    // (none)
    // Amplitude(s) for diagram number 13550
    FFV1_0( w_fp[450], w_fp[2], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 13551 OF 15495 ***
    // Wavefunction(s) for diagram number 13551
    // (none)
    // Amplitude(s) for diagram number 13551
    VVVV1_0( w_fp[0], w_fp[1], w_fp[154], w_fp[88], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[154], w_fp[88], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[0], w_fp[1], w_fp[154], w_fp[88], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13552 OF 15495 ***
    // Wavefunction(s) for diagram number 13552
    // (none)
    // Amplitude(s) for diagram number 13552
    VVV1_0( w_fp[1], w_fp[88], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13553 OF 15495 ***
    // Wavefunction(s) for diagram number 13553
    // (none)
    // Amplitude(s) for diagram number 13553
    VVV1_0( w_fp[1], w_fp[154], w_fp[609], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13554 OF 15495 ***
    // Wavefunction(s) for diagram number 13554
    // (none)
    // Amplitude(s) for diagram number 13554
    FFV1_0( w_fp[171], w_fp[92], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13555 OF 15495 ***
    // Wavefunction(s) for diagram number 13555
    // (none)
    // Amplitude(s) for diagram number 13555
    FFV1_0( w_fp[450], w_fp[118], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13556 OF 15495 ***
    // Wavefunction(s) for diagram number 13556
    // (none)
    // Amplitude(s) for diagram number 13556
    VVV1_0( w_fp[22], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVV1_0( w_fp[21], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0( w_fp[20], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 13557 OF 15495 ***
    // Wavefunction(s) for diagram number 13557
    // (none)
    // Amplitude(s) for diagram number 13557
    FFV1_0( w_fp[171], w_fp[2], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13558 OF 15495 ***
    // Wavefunction(s) for diagram number 13558
    // (none)
    // Amplitude(s) for diagram number 13558
    FFV1_0( w_fp[256], w_fp[2], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13559 OF 15495 ***
    // Wavefunction(s) for diagram number 13559
    // (none)
    // Amplitude(s) for diagram number 13559
    VVV1_0( w_fp[550], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[480], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    VVV1_0( w_fp[558], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 13560 OF 15495 ***
    // Wavefunction(s) for diagram number 13560
    // (none)
    // Amplitude(s) for diagram number 13560
    VVV1_0( w_fp[0], w_fp[367], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[368], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[369], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13561 OF 15495 ***
    // Wavefunction(s) for diagram number 13561
    // (none)
    // Amplitude(s) for diagram number 13561
    VVV1_0( w_fp[363], w_fp[6], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13562 OF 15495 ***
    // Wavefunction(s) for diagram number 13562
    // (none)
    // Amplitude(s) for diagram number 13562
    FFV1_0( w_fp[538], w_fp[2], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];

    // *** DIAGRAM 13563 OF 15495 ***
    // Wavefunction(s) for diagram number 13563
    // (none)
    // Amplitude(s) for diagram number 13563
    FFV1_0( w_fp[660], w_fp[98], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13564 OF 15495 ***
    // Wavefunction(s) for diagram number 13564
    // (none)
    // Amplitude(s) for diagram number 13564
    FFV1_0( w_fp[538], w_fp[98], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13565 OF 15495 ***
    // Wavefunction(s) for diagram number 13565
    // (none)
    // Amplitude(s) for diagram number 13565
    FFV1_0( w_fp[660], w_fp[2], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[356] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];

    // *** DIAGRAM 13566 OF 15495 ***
    // Wavefunction(s) for diagram number 13566
    // (none)
    // Amplitude(s) for diagram number 13566
    VVV1_0( w_fp[1], w_fp[104], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13567 OF 15495 ***
    // Wavefunction(s) for diagram number 13567
    // (none)
    // Amplitude(s) for diagram number 13567
    FFV1_0( w_fp[527], w_fp[2], w_fp[371], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[372], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[373], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13568 OF 15495 ***
    // Wavefunction(s) for diagram number 13568
    // (none)
    // Amplitude(s) for diagram number 13568
    FFV1_0( w_fp[256], w_fp[245], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 13569 OF 15495 ***
    // Wavefunction(s) for diagram number 13569
    // (none)
    // Amplitude(s) for diagram number 13569
    FFV1_0( w_fp[256], w_fp[2], w_fp[613], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13570 OF 15495 ***
    // Wavefunction(s) for diagram number 13570
    // (none)
    // Amplitude(s) for diagram number 13570
    VVVV1_0( w_fp[518], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    VVVV3_0( w_fp[518], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    VVVV4_0( w_fp[518], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 13571 OF 15495 ***
    // Wavefunction(s) for diagram number 13571
    // (none)
    // Amplitude(s) for diagram number 13571
    VVV1_0( w_fp[154], w_fp[6], w_fp[611], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 13572 OF 15495 ***
    // Wavefunction(s) for diagram number 13572
    // (none)
    // Amplitude(s) for diagram number 13572
    VVV1_0( w_fp[1], w_fp[154], w_fp[613], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 13573 OF 15495 ***
    // Wavefunction(s) for diagram number 13573
    // (none)
    // Amplitude(s) for diagram number 13573
    FFV1_0( w_fp[170], w_fp[2], w_fp[611], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13574 OF 15495 ***
    // Wavefunction(s) for diagram number 13574
    // (none)
    // Amplitude(s) for diagram number 13574
    FFV1_0( w_fp[170], w_fp[245], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 13575 OF 15495 ***
    // Wavefunction(s) for diagram number 13575
    // (none)
    // Amplitude(s) for diagram number 13575
    FFV1_0( w_fp[658], w_fp[98], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13576 OF 15495 ***
    // Wavefunction(s) for diagram number 13576
    // (none)
    // Amplitude(s) for diagram number 13576
    FFV1_0( w_fp[256], w_fp[110], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13577 OF 15495 ***
    // Wavefunction(s) for diagram number 13577
    // (none)
    // Amplitude(s) for diagram number 13577
    FFV1_0( w_fp[658], w_fp[2], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];

    // *** DIAGRAM 13578 OF 15495 ***
    // Wavefunction(s) for diagram number 13578
    // (none)
    // Amplitude(s) for diagram number 13578
    FFV1_0( w_fp[256], w_fp[2], w_fp[616], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13579 OF 15495 ***
    // Wavefunction(s) for diagram number 13579
    // (none)
    // Amplitude(s) for diagram number 13579
    VVVV1_0( w_fp[0], w_fp[363], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[363], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[363], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13580 OF 15495 ***
    // Wavefunction(s) for diagram number 13580
    // (none)
    // Amplitude(s) for diagram number 13580
    VVV1_0( w_fp[154], w_fp[6], w_fp[61], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 13581 OF 15495 ***
    // Wavefunction(s) for diagram number 13581
    // (none)
    // Amplitude(s) for diagram number 13581
    VVV1_0( w_fp[363], w_fp[6], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13582 OF 15495 ***
    // Wavefunction(s) for diagram number 13582
    // (none)
    // Amplitude(s) for diagram number 13582
    FFV1_0( w_fp[170], w_fp[2], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13583 OF 15495 ***
    // Wavefunction(s) for diagram number 13583
    // (none)
    // Amplitude(s) for diagram number 13583
    FFV1_0( w_fp[535], w_fp[2], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 13584 OF 15495 ***
    // Wavefunction(s) for diagram number 13584
    // (none)
    // Amplitude(s) for diagram number 13584
    VVVV1_0( w_fp[0], w_fp[1], w_fp[154], w_fp[104], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[154], w_fp[104], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[0], w_fp[1], w_fp[154], w_fp[104], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13585 OF 15495 ***
    // Wavefunction(s) for diagram number 13585
    // (none)
    // Amplitude(s) for diagram number 13585
    VVV1_0( w_fp[1], w_fp[104], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13586 OF 15495 ***
    // Wavefunction(s) for diagram number 13586
    // (none)
    // Amplitude(s) for diagram number 13586
    VVV1_0( w_fp[1], w_fp[154], w_fp[616], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13587 OF 15495 ***
    // Wavefunction(s) for diagram number 13587
    // (none)
    // Amplitude(s) for diagram number 13587
    FFV1_0( w_fp[170], w_fp[110], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13588 OF 15495 ***
    // Wavefunction(s) for diagram number 13588
    // (none)
    // Amplitude(s) for diagram number 13588
    FFV1_0( w_fp[535], w_fp[98], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13589 OF 15495 ***
    // Wavefunction(s) for diagram number 13589
    // (none)
    // Amplitude(s) for diagram number 13589
    VVV1_0( w_fp[19], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVV1_0( w_fp[18], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    VVV1_0( w_fp[17], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 13590 OF 15495 ***
    // Wavefunction(s) for diagram number 13590
    // (none)
    // Amplitude(s) for diagram number 13590
    FFV1_0( w_fp[170], w_fp[2], w_fp[19], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[18], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[17], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13591 OF 15495 ***
    // Wavefunction(s) for diagram number 13591
    // (none)
    // Amplitude(s) for diagram number 13591
    FFV1_0( w_fp[256], w_fp[2], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[529], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13592 OF 15495 ***
    // Wavefunction(s) for diagram number 13592
    // (none)
    // Amplitude(s) for diagram number 13592
    VVV1_0( w_fp[486], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[547], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    VVV1_0( w_fp[529], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 13593 OF 15495 ***
    // Wavefunction(s) for diagram number 13593
    // (none)
    // Amplitude(s) for diagram number 13593
    VVV1_0( w_fp[0], w_fp[371], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[372], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[373], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 13594 OF 15495 ***
    // Wavefunction(s) for diagram number 13594
    // (none)
    // Amplitude(s) for diagram number 13594
    VVV1_0( w_fp[361], w_fp[4], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13595 OF 15495 ***
    // Wavefunction(s) for diagram number 13595
    // (none)
    // Amplitude(s) for diagram number 13595
    FFV1_0( w_fp[509], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];

    // *** DIAGRAM 13596 OF 15495 ***
    // Wavefunction(s) for diagram number 13596
    // (none)
    // Amplitude(s) for diagram number 13596
    FFV1_0( w_fp[660], w_fp[128], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13597 OF 15495 ***
    // Wavefunction(s) for diagram number 13597
    // (none)
    // Amplitude(s) for diagram number 13597
    FFV1_0( w_fp[509], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13598 OF 15495 ***
    // Wavefunction(s) for diagram number 13598
    // (none)
    // Amplitude(s) for diagram number 13598
    FFV1_0( w_fp[660], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[332] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 13599 OF 15495 ***
    // Wavefunction(s) for diagram number 13599
    // (none)
    // Amplitude(s) for diagram number 13599
    VVV1_0( w_fp[1], w_fp[130], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 13600 OF 15495 ***
    // Wavefunction(s) for diagram number 13600
    // (none)
    // Amplitude(s) for diagram number 13600
    FFV1_0( w_fp[527], w_fp[2], w_fp[291], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[527], w_fp[2], w_fp[248], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[264], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

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
