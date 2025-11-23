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
  diagramgroup97( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 121 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 138 );
    retrieveWf( wfs, w_cx, nevt, 155 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 236 );
    retrieveWf( wfs, w_cx, nevt, 237 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 465 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 517 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 587 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 600 );
    retrieveWf( wfs, w_cx, nevt, 618 );
    retrieveWf( wfs, w_cx, nevt, 619 );
    retrieveWf( wfs, w_cx, nevt, 621 );
    retrieveWf( wfs, w_cx, nevt, 622 );
    retrieveWf( wfs, w_cx, nevt, 625 );
    retrieveWf( wfs, w_cx, nevt, 626 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 663 );
    retrieveWf( wfs, w_cx, nevt, 664 );
    retrieveWf( wfs, w_cx, nevt, 665 );
    retrieveWf( wfs, w_cx, nevt, 673 );
    retrieveWf( wfs, w_cx, nevt, 675 );
    retrieveWf( wfs, w_cx, nevt, 680 );
    retrieveWf( wfs, w_cx, nevt, 697 );
#endif
#endif

    // *** DIAGRAM 9601 OF 15495 ***
    // Wavefunction(s) for diagram number 9601
    // (none)
    // Amplitude(s) for diagram number 9601
    FFV1_0( w_fp[662], w_fp[2], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9602 OF 15495 ***
    // Wavefunction(s) for diagram number 9602
    // (none)
    // Amplitude(s) for diagram number 9602
    VVV1_0( w_fp[518], w_fp[534], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];

    // *** DIAGRAM 9603 OF 15495 ***
    // Wavefunction(s) for diagram number 9603
    // (none)
    // Amplitude(s) for diagram number 9603
    FFV1_0( w_fp[121], w_fp[2], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9604 OF 15495 ***
    // Wavefunction(s) for diagram number 9604
    // (none)
    // Amplitude(s) for diagram number 9604
    FFV1_0( w_fp[157], w_fp[244], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9605 OF 15495 ***
    // Wavefunction(s) for diagram number 9605
    // (none)
    // Amplitude(s) for diagram number 9605
    VVV1_0( w_fp[559], w_fp[534], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];

    // *** DIAGRAM 9606 OF 15495 ***
    // Wavefunction(s) for diagram number 9606
    // (none)
    // Amplitude(s) for diagram number 9606
    FFV1_0( w_fp[120], w_fp[2], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9607 OF 15495 ***
    // Wavefunction(s) for diagram number 9607
    // (none)
    // Amplitude(s) for diagram number 9607
    FFV1_0( w_fp[157], w_fp[98], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9608 OF 15495 ***
    // Wavefunction(s) for diagram number 9608
    // (none)
    // Amplitude(s) for diagram number 9608
    VVV1_0( w_fp[0], w_fp[534], w_fp[59], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];

    // *** DIAGRAM 9609 OF 15495 ***
    // Wavefunction(s) for diagram number 9609
    // (none)
    // Amplitude(s) for diagram number 9609
    FFV1_0( w_fp[120], w_fp[244], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];

    // *** DIAGRAM 9610 OF 15495 ***
    // Wavefunction(s) for diagram number 9610
    // (none)
    // Amplitude(s) for diagram number 9610
    FFV1_0( w_fp[121], w_fp[98], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 9611 OF 15495 ***
    // Wavefunction(s) for diagram number 9611
    // (none)
    // Amplitude(s) for diagram number 9611
    FFV1_0( w_fp[157], w_fp[2], w_fp[465], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[463], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[618], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];

    // *** DIAGRAM 9612 OF 15495 ***
    // Wavefunction(s) for diagram number 9612
    // (none)
    // Amplitude(s) for diagram number 9612
    FFV1_0( w_fp[497], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9613 OF 15495 ***
    // Wavefunction(s) for diagram number 9613
    // (none)
    // Amplitude(s) for diagram number 9613
    FFV1_0( w_fp[664], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9614 OF 15495 ***
    // Wavefunction(s) for diagram number 9614
    // (none)
    // Amplitude(s) for diagram number 9614
    VVV1_0( w_fp[115], w_fp[7], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9615 OF 15495 ***
    // Wavefunction(s) for diagram number 9615
    // (none)
    // Amplitude(s) for diagram number 9615
    FFV1_0( w_fp[664], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];

    // *** DIAGRAM 9616 OF 15495 ***
    // Wavefunction(s) for diagram number 9616
    // (none)
    // Amplitude(s) for diagram number 9616
    VVV1_0( w_fp[4], w_fp[116], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9617 OF 15495 ***
    // Wavefunction(s) for diagram number 9617
    // (none)
    // Amplitude(s) for diagram number 9617
    FFV1_0( w_fp[497], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9618 OF 15495 ***
    // Wavefunction(s) for diagram number 9618
    // (none)
    // Amplitude(s) for diagram number 9618
    FFV1_0( w_fp[662], w_fp[2], w_fp[237], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[662], w_fp[2], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[662], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9619 OF 15495 ***
    // Wavefunction(s) for diagram number 9619
    // (none)
    // Amplitude(s) for diagram number 9619
    VVVV1_0( w_fp[559], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    VVVV3_0( w_fp[559], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVVV4_0( w_fp[559], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9620 OF 15495 ***
    // Wavefunction(s) for diagram number 9620
    // (none)
    // Amplitude(s) for diagram number 9620
    VVV1_0( w_fp[534], w_fp[7], w_fp[585], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9621 OF 15495 ***
    // Wavefunction(s) for diagram number 9621
    // (none)
    // Amplitude(s) for diagram number 9621
    VVV1_0( w_fp[534], w_fp[4], w_fp[619], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9622 OF 15495 ***
    // Wavefunction(s) for diagram number 9622
    FFV1_1( w_fp[2], w_fp[559], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[664] );
    // Amplitude(s) for diagram number 9622
    FFV1_0( w_fp[94], w_fp[664], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];

    // *** DIAGRAM 9623 OF 15495 ***
    // Wavefunction(s) for diagram number 9623
    // (none)
    // Amplitude(s) for diagram number 9623
    FFV1_0( w_fp[94], w_fp[2], w_fp[619], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9624 OF 15495 ***
    // Wavefunction(s) for diagram number 9624
    // (none)
    // Amplitude(s) for diagram number 9624
    FFV1_0( w_fp[77], w_fp[664], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];

    // *** DIAGRAM 9625 OF 15495 ***
    // Wavefunction(s) for diagram number 9625
    // (none)
    // Amplitude(s) for diagram number 9625
    FFV1_0( w_fp[77], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9626 OF 15495 ***
    // Wavefunction(s) for diagram number 9626
    // (none)
    // Amplitude(s) for diagram number 9626
    VVVV1_0( w_fp[0], w_fp[534], w_fp[115], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[534], w_fp[115], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[534], w_fp[115], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];

    // *** DIAGRAM 9627 OF 15495 ***
    // Wavefunction(s) for diagram number 9627
    // (none)
    // Amplitude(s) for diagram number 9627
    VVV1_0( w_fp[115], w_fp[7], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9628 OF 15495 ***
    // Wavefunction(s) for diagram number 9628
    // (none)
    // Amplitude(s) for diagram number 9628
    VVV1_0( w_fp[534], w_fp[7], w_fp[621], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9629 OF 15495 ***
    // Wavefunction(s) for diagram number 9629
    // (none)
    // Amplitude(s) for diagram number 9629
    VVVV1_0( w_fp[0], w_fp[534], w_fp[4], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[534], w_fp[4], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[534], w_fp[4], w_fp[116], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9630 OF 15495 ***
    // Wavefunction(s) for diagram number 9630
    // (none)
    // Amplitude(s) for diagram number 9630
    VVV1_0( w_fp[4], w_fp[116], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9631 OF 15495 ***
    // Wavefunction(s) for diagram number 9631
    // (none)
    // Amplitude(s) for diagram number 9631
    VVV1_0( w_fp[534], w_fp[4], w_fp[622], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9632 OF 15495 ***
    // Wavefunction(s) for diagram number 9632
    // (none)
    // Amplitude(s) for diagram number 9632
    VVV1_0( w_fp[0], w_fp[534], w_fp[237], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[534], w_fp[8], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[534], w_fp[114], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 9633 OF 15495 ***
    // Wavefunction(s) for diagram number 9633
    // (none)
    // Amplitude(s) for diagram number 9633
    FFV1_0( w_fp[680], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9634 OF 15495 ***
    // Wavefunction(s) for diagram number 9634
    FFV1_1( w_fp[244], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[121] );
    // Amplitude(s) for diagram number 9634
    FFV1_0( w_fp[94], w_fp[121], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9635 OF 15495 ***
    // Wavefunction(s) for diagram number 9635
    // (none)
    // Amplitude(s) for diagram number 9635
    FFV1_0( w_fp[680], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 9636 OF 15495 ***
    // Wavefunction(s) for diagram number 9636
    // (none)
    // Amplitude(s) for diagram number 9636
    FFV1_0( w_fp[94], w_fp[2], w_fp[622], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9637 OF 15495 ***
    // Wavefunction(s) for diagram number 9637
    // (none)
    // Amplitude(s) for diagram number 9637
    FFV1_0( w_fp[675], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9638 OF 15495 ***
    // Wavefunction(s) for diagram number 9638
    // (none)
    // Amplitude(s) for diagram number 9638
    FFV1_0( w_fp[77], w_fp[121], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9639 OF 15495 ***
    // Wavefunction(s) for diagram number 9639
    // (none)
    // Amplitude(s) for diagram number 9639
    FFV1_0( w_fp[675], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];

    // *** DIAGRAM 9640 OF 15495 ***
    // Wavefunction(s) for diagram number 9640
    // (none)
    // Amplitude(s) for diagram number 9640
    FFV1_0( w_fp[77], w_fp[2], w_fp[621], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9641 OF 15495 ***
    // Wavefunction(s) for diagram number 9641
    // (none)
    // Amplitude(s) for diagram number 9641
    VVV1_0( w_fp[517], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[472], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    VVV1_0( w_fp[445], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 9642 OF 15495 ***
    // Wavefunction(s) for diagram number 9642
    // (none)
    // Amplitude(s) for diagram number 9642
    FFV1_0( w_fp[77], w_fp[2], w_fp[517], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[472], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[445], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9643 OF 15495 ***
    // Wavefunction(s) for diagram number 9643
    // (none)
    // Amplitude(s) for diagram number 9643
    VVV1_0( w_fp[479], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVV1_0( w_fp[474], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    VVV1_0( w_fp[444], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 9644 OF 15495 ***
    // Wavefunction(s) for diagram number 9644
    // (none)
    // Amplitude(s) for diagram number 9644
    FFV1_0( w_fp[94], w_fp[2], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9645 OF 15495 ***
    // Wavefunction(s) for diagram number 9645
    // (none)
    // Amplitude(s) for diagram number 9645
    FFV1_0( w_fp[497], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9646 OF 15495 ***
    // Wavefunction(s) for diagram number 9646
    // (none)
    // Amplitude(s) for diagram number 9646
    FFV1_0( w_fp[663], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9647 OF 15495 ***
    // Wavefunction(s) for diagram number 9647
    // (none)
    // Amplitude(s) for diagram number 9647
    VVV1_0( w_fp[124], w_fp[6], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9648 OF 15495 ***
    // Wavefunction(s) for diagram number 9648
    // (none)
    // Amplitude(s) for diagram number 9648
    FFV1_0( w_fp[663], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];

    // *** DIAGRAM 9649 OF 15495 ***
    // Wavefunction(s) for diagram number 9649
    // (none)
    // Amplitude(s) for diagram number 9649
    VVV1_0( w_fp[4], w_fp[125], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9650 OF 15495 ***
    // Wavefunction(s) for diagram number 9650
    // (none)
    // Amplitude(s) for diagram number 9650
    FFV1_0( w_fp[497], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[478] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];

    // *** DIAGRAM 9651 OF 15495 ***
    // Wavefunction(s) for diagram number 9651
    // (none)
    // Amplitude(s) for diagram number 9651
    FFV1_0( w_fp[662], w_fp[2], w_fp[236], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[662], w_fp[2], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[662], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9652 OF 15495 ***
    // Wavefunction(s) for diagram number 9652
    // (none)
    // Amplitude(s) for diagram number 9652
    VVVV1_0( w_fp[522], w_fp[534], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    VVVV3_0( w_fp[522], w_fp[534], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    VVVV4_0( w_fp[522], w_fp[534], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 9653 OF 15495 ***
    // Wavefunction(s) for diagram number 9653
    // (none)
    // Amplitude(s) for diagram number 9653
    VVV1_0( w_fp[534], w_fp[6], w_fp[453], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 9654 OF 15495 ***
    // Wavefunction(s) for diagram number 9654
    // (none)
    // Amplitude(s) for diagram number 9654
    VVV1_0( w_fp[534], w_fp[4], w_fp[505], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 9655 OF 15495 ***
    // Wavefunction(s) for diagram number 9655
    FFV1_1( w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[663] );
    // Amplitude(s) for diagram number 9655
    FFV1_0( w_fp[94], w_fp[663], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 9656 OF 15495 ***
    // Wavefunction(s) for diagram number 9656
    // (none)
    // Amplitude(s) for diagram number 9656
    FFV1_0( w_fp[94], w_fp[2], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9657 OF 15495 ***
    // Wavefunction(s) for diagram number 9657
    // (none)
    // Amplitude(s) for diagram number 9657
    FFV1_0( w_fp[39], w_fp[663], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 9658 OF 15495 ***
    // Wavefunction(s) for diagram number 9658
    // (none)
    // Amplitude(s) for diagram number 9658
    FFV1_0( w_fp[39], w_fp[2], w_fp[453], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9659 OF 15495 ***
    // Wavefunction(s) for diagram number 9659
    // (none)
    // Amplitude(s) for diagram number 9659
    VVVV1_0( w_fp[0], w_fp[534], w_fp[124], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[534], w_fp[124], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[534], w_fp[124], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];

    // *** DIAGRAM 9660 OF 15495 ***
    // Wavefunction(s) for diagram number 9660
    // (none)
    // Amplitude(s) for diagram number 9660
    VVV1_0( w_fp[124], w_fp[6], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];

    // *** DIAGRAM 9661 OF 15495 ***
    // Wavefunction(s) for diagram number 9661
    // (none)
    // Amplitude(s) for diagram number 9661
    VVV1_0( w_fp[534], w_fp[6], w_fp[625], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];

    // *** DIAGRAM 9662 OF 15495 ***
    // Wavefunction(s) for diagram number 9662
    // (none)
    // Amplitude(s) for diagram number 9662
    VVVV1_0( w_fp[0], w_fp[534], w_fp[4], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[534], w_fp[4], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[534], w_fp[4], w_fp[125], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];

    // *** DIAGRAM 9663 OF 15495 ***
    // Wavefunction(s) for diagram number 9663
    // (none)
    // Amplitude(s) for diagram number 9663
    VVV1_0( w_fp[4], w_fp[125], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];

    // *** DIAGRAM 9664 OF 15495 ***
    // Wavefunction(s) for diagram number 9664
    // (none)
    // Amplitude(s) for diagram number 9664
    VVV1_0( w_fp[534], w_fp[4], w_fp[626], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];

    // *** DIAGRAM 9665 OF 15495 ***
    // Wavefunction(s) for diagram number 9665
    // (none)
    // Amplitude(s) for diagram number 9665
    VVV1_0( w_fp[0], w_fp[534], w_fp[236], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[534], w_fp[56], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[534], w_fp[55], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];

    // *** DIAGRAM 9666 OF 15495 ***
    // Wavefunction(s) for diagram number 9666
    // (none)
    // Amplitude(s) for diagram number 9666
    FFV1_0( w_fp[680], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9667 OF 15495 ***
    // Wavefunction(s) for diagram number 9667
    FFV1_1( w_fp[122], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[675] );
    // Amplitude(s) for diagram number 9667
    FFV1_0( w_fp[94], w_fp[675], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9668 OF 15495 ***
    // Wavefunction(s) for diagram number 9668
    // (none)
    // Amplitude(s) for diagram number 9668
    FFV1_0( w_fp[680], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[475] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];

    // *** DIAGRAM 9669 OF 15495 ***
    // Wavefunction(s) for diagram number 9669
    // (none)
    // Amplitude(s) for diagram number 9669
    FFV1_0( w_fp[94], w_fp[2], w_fp[626], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9670 OF 15495 ***
    // Wavefunction(s) for diagram number 9670
    // (none)
    // Amplitude(s) for diagram number 9670
    FFV1_0( w_fp[673], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9671 OF 15495 ***
    // Wavefunction(s) for diagram number 9671
    // (none)
    // Amplitude(s) for diagram number 9671
    FFV1_0( w_fp[39], w_fp[675], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9672 OF 15495 ***
    // Wavefunction(s) for diagram number 9672
    // (none)
    // Amplitude(s) for diagram number 9672
    FFV1_0( w_fp[673], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];

    // *** DIAGRAM 9673 OF 15495 ***
    // Wavefunction(s) for diagram number 9673
    // (none)
    // Amplitude(s) for diagram number 9673
    FFV1_0( w_fp[39], w_fp[2], w_fp[625], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9674 OF 15495 ***
    // Wavefunction(s) for diagram number 9674
    // (none)
    // Amplitude(s) for diagram number 9674
    VVV1_0( w_fp[475], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    VVV1_0( w_fp[553], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    VVV1_0( w_fp[592], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];

    // *** DIAGRAM 9675 OF 15495 ***
    // Wavefunction(s) for diagram number 9675
    // (none)
    // Amplitude(s) for diagram number 9675
    FFV1_0( w_fp[39], w_fp[2], w_fp[475], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[39], w_fp[2], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[39], w_fp[2], w_fp[592], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9676 OF 15495 ***
    // Wavefunction(s) for diagram number 9676
    // (none)
    // Amplitude(s) for diagram number 9676
    VVV1_0( w_fp[523], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    VVV1_0( w_fp[580], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    VVV1_0( w_fp[528], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];

    // *** DIAGRAM 9677 OF 15495 ***
    // Wavefunction(s) for diagram number 9677
    // (none)
    // Amplitude(s) for diagram number 9677
    FFV1_0( w_fp[94], w_fp[2], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9678 OF 15495 ***
    // Wavefunction(s) for diagram number 9678
    // (none)
    // Amplitude(s) for diagram number 9678
    FFV1_0( w_fp[497], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9679 OF 15495 ***
    // Wavefunction(s) for diagram number 9679
    // (none)
    // Amplitude(s) for diagram number 9679
    FFV1_0( w_fp[665], w_fp[128], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9680 OF 15495 ***
    // Wavefunction(s) for diagram number 9680
    // (none)
    // Amplitude(s) for diagram number 9680
    VVV1_0( w_fp[130], w_fp[5], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9681 OF 15495 ***
    // Wavefunction(s) for diagram number 9681
    // (none)
    // Amplitude(s) for diagram number 9681
    FFV1_0( w_fp[665], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 9682 OF 15495 ***
    // Wavefunction(s) for diagram number 9682
    // (none)
    // Amplitude(s) for diagram number 9682
    VVV1_0( w_fp[4], w_fp[131], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9683 OF 15495 ***
    // Wavefunction(s) for diagram number 9683
    // (none)
    // Amplitude(s) for diagram number 9683
    FFV1_0( w_fp[497], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9684 OF 15495 ***
    // Wavefunction(s) for diagram number 9684
    // (none)
    // Amplitude(s) for diagram number 9684
    FFV1_0( w_fp[662], w_fp[2], w_fp[155], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[662], w_fp[2], w_fp[138], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[662], w_fp[2], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9685 OF 15495 ***
    // Wavefunction(s) for diagram number 9685
    // (none)
    // Amplitude(s) for diagram number 9685
    VVVV1_0( w_fp[516], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVVV3_0( w_fp[516], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    VVVV4_0( w_fp[516], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 9686 OF 15495 ***
    // Wavefunction(s) for diagram number 9686
    // (none)
    // Amplitude(s) for diagram number 9686
    VVV1_0( w_fp[534], w_fp[5], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 9687 OF 15495 ***
    // Wavefunction(s) for diagram number 9687
    // (none)
    // Amplitude(s) for diagram number 9687
    VVV1_0( w_fp[534], w_fp[4], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 9688 OF 15495 ***
    // Wavefunction(s) for diagram number 9688
    FFV1_1( w_fp[2], w_fp[516], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[497] );
    // Amplitude(s) for diagram number 9688
    FFV1_0( w_fp[94], w_fp[497], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 9689 OF 15495 ***
    // Wavefunction(s) for diagram number 9689
    // (none)
    // Amplitude(s) for diagram number 9689
    FFV1_0( w_fp[94], w_fp[2], w_fp[587], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9690 OF 15495 ***
    // Wavefunction(s) for diagram number 9690
    // (none)
    // Amplitude(s) for diagram number 9690
    FFV1_0( w_fp[29], w_fp[497], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 9691 OF 15495 ***
    // Wavefunction(s) for diagram number 9691
    // (none)
    // Amplitude(s) for diagram number 9691
    FFV1_0( w_fp[29], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9692 OF 15495 ***
    // Wavefunction(s) for diagram number 9692
    // (none)
    // Amplitude(s) for diagram number 9692
    VVVV1_0( w_fp[0], w_fp[534], w_fp[130], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[534], w_fp[130], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[534], w_fp[130], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 9693 OF 15495 ***
    // Wavefunction(s) for diagram number 9693
    // (none)
    // Amplitude(s) for diagram number 9693
    VVV1_0( w_fp[130], w_fp[5], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 9694 OF 15495 ***
    // Wavefunction(s) for diagram number 9694
    // (none)
    // Amplitude(s) for diagram number 9694
    VVV1_0( w_fp[534], w_fp[5], w_fp[525], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 9695 OF 15495 ***
    // Wavefunction(s) for diagram number 9695
    // (none)
    // Amplitude(s) for diagram number 9695
    VVVV1_0( w_fp[0], w_fp[534], w_fp[4], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[534], w_fp[4], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[534], w_fp[4], w_fp[131], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9696 OF 15495 ***
    // Wavefunction(s) for diagram number 9696
    // (none)
    // Amplitude(s) for diagram number 9696
    VVV1_0( w_fp[4], w_fp[131], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 9697 OF 15495 ***
    // Wavefunction(s) for diagram number 9697
    // (none)
    // Amplitude(s) for diagram number 9697
    VVV1_0( w_fp[534], w_fp[4], w_fp[600], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 9698 OF 15495 ***
    // Wavefunction(s) for diagram number 9698
    // (none)
    // Amplitude(s) for diagram number 9698
    VVV1_0( w_fp[0], w_fp[534], w_fp[155], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[534], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[534], w_fp[238], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 9699 OF 15495 ***
    // Wavefunction(s) for diagram number 9699
    // (none)
    // Amplitude(s) for diagram number 9699
    FFV1_0( w_fp[680], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9700 OF 15495 ***
    // Wavefunction(s) for diagram number 9700
    FFV1_1( w_fp[128], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[697] );
    // Amplitude(s) for diagram number 9700
    FFV1_0( w_fp[94], w_fp[697], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 121 );
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 663 );
    storeWf( wfs, w_cx, nevt, 664 );
    storeWf( wfs, w_cx, nevt, 675 );
    storeWf( wfs, w_cx, nevt, 697 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
