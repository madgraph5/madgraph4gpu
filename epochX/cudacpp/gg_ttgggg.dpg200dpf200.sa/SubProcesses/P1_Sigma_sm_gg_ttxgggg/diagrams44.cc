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
  diagramgroup44( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 48 );
    retrieveWf( wfs, w_cx, nevt, 49 );
    retrieveWf( wfs, w_cx, nevt, 50 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 93 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 112 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 142 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 172 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 195 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 262 );
    retrieveWf( wfs, w_cx, nevt, 276 );
    retrieveWf( wfs, w_cx, nevt, 303 );
    retrieveWf( wfs, w_cx, nevt, 353 );
    retrieveWf( wfs, w_cx, nevt, 356 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 376 );
    retrieveWf( wfs, w_cx, nevt, 377 );
    retrieveWf( wfs, w_cx, nevt, 378 );
    retrieveWf( wfs, w_cx, nevt, 379 );
    retrieveWf( wfs, w_cx, nevt, 380 );
    retrieveWf( wfs, w_cx, nevt, 381 );
    retrieveWf( wfs, w_cx, nevt, 382 );
    retrieveWf( wfs, w_cx, nevt, 386 );
    retrieveWf( wfs, w_cx, nevt, 387 );
    retrieveWf( wfs, w_cx, nevt, 388 );
    retrieveWf( wfs, w_cx, nevt, 400 );
    retrieveWf( wfs, w_cx, nevt, 401 );
    retrieveWf( wfs, w_cx, nevt, 402 );
    retrieveWf( wfs, w_cx, nevt, 403 );
    retrieveWf( wfs, w_cx, nevt, 404 );
    retrieveWf( wfs, w_cx, nevt, 408 );
    retrieveWf( wfs, w_cx, nevt, 409 );
    retrieveWf( wfs, w_cx, nevt, 410 );
    retrieveWf( wfs, w_cx, nevt, 432 );
    retrieveWf( wfs, w_cx, nevt, 433 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 439 );
    retrieveWf( wfs, w_cx, nevt, 440 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 446 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 460 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 464 );
    retrieveWf( wfs, w_cx, nevt, 465 );
    retrieveWf( wfs, w_cx, nevt, 466 );
    retrieveWf( wfs, w_cx, nevt, 467 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 482 );
    retrieveWf( wfs, w_cx, nevt, 484 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 513 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 541 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
#endif
#endif

    // *** DIAGRAM 8601 OF 15495 ***
    // Wavefunction(s) for diagram number 8601
    // (none)
    // Amplitude(s) for diagram number 8601
    FFV1_0( w_fp[184], w_fp[437], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8602 OF 15495 ***
    // Wavefunction(s) for diagram number 8602
    // (none)
    // Amplitude(s) for diagram number 8602
    FFV1_0( w_fp[184], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];

    // *** DIAGRAM 8603 OF 15495 ***
    // Wavefunction(s) for diagram number 8603
    // (none)
    // Amplitude(s) for diagram number 8603
    FFV1_0( w_fp[3], w_fp[437], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 8604 OF 15495 ***
    // Wavefunction(s) for diagram number 8604
    // (none)
    // Amplitude(s) for diagram number 8604
    FFV1_0( w_fp[488], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];

    // *** DIAGRAM 8605 OF 15495 ***
    // Wavefunction(s) for diagram number 8605
    // (none)
    // Amplitude(s) for diagram number 8605
    VVV1_0( w_fp[359], w_fp[4], w_fp[451], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 8606 OF 15495 ***
    // Wavefunction(s) for diagram number 8606
    // (none)
    // Amplitude(s) for diagram number 8606
    FFV1_0( w_fp[3], w_fp[452], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8607 OF 15495 ***
    // Wavefunction(s) for diagram number 8607
    // (none)
    // Amplitude(s) for diagram number 8607
    FFV1_0( w_fp[184], w_fp[571], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];

    // *** DIAGRAM 8608 OF 15495 ***
    // Wavefunction(s) for diagram number 8608
    // (none)
    // Amplitude(s) for diagram number 8608
    FFV1_0( w_fp[184], w_fp[452], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];

    // *** DIAGRAM 8609 OF 15495 ***
    // Wavefunction(s) for diagram number 8609
    // (none)
    // Amplitude(s) for diagram number 8609
    FFV1_0( w_fp[3], w_fp[571], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8610 OF 15495 ***
    // Wavefunction(s) for diagram number 8610
    // (none)
    // Amplitude(s) for diagram number 8610
    VVV1_0( w_fp[1], w_fp[115], w_fp[451], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];

    // *** DIAGRAM 8611 OF 15495 ***
    // Wavefunction(s) for diagram number 8611
    // (none)
    // Amplitude(s) for diagram number 8611
    FFV1_0( w_fp[3], w_fp[453], w_fp[262], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[353], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[72], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];

    // *** DIAGRAM 8612 OF 15495 ***
    // Wavefunction(s) for diagram number 8612
    // (none)
    // Amplitude(s) for diagram number 8612
    VVV1_0( w_fp[359], w_fp[4], w_fp[580], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];

    // *** DIAGRAM 8613 OF 15495 ***
    // Wavefunction(s) for diagram number 8613
    // (none)
    // Amplitude(s) for diagram number 8613
    FFV1_0( w_fp[516], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8614 OF 15495 ***
    // Wavefunction(s) for diagram number 8614
    // (none)
    // Amplitude(s) for diagram number 8614
    FFV1_0( w_fp[528], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 8615 OF 15495 ***
    // Wavefunction(s) for diagram number 8615
    // (none)
    // Amplitude(s) for diagram number 8615
    FFV1_0( w_fp[516], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];

    // *** DIAGRAM 8616 OF 15495 ***
    // Wavefunction(s) for diagram number 8616
    // (none)
    // Amplitude(s) for diagram number 8616
    FFV1_0( w_fp[528], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8617 OF 15495 ***
    // Wavefunction(s) for diagram number 8617
    // (none)
    // Amplitude(s) for diagram number 8617
    VVV1_0( w_fp[1], w_fp[115], w_fp[580], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];

    // *** DIAGRAM 8618 OF 15495 ***
    // Wavefunction(s) for diagram number 8618
    // (none)
    // Amplitude(s) for diagram number 8618
    FFV1_0( w_fp[449], w_fp[2], w_fp[262], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[353], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[72], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 8619 OF 15495 ***
    // Wavefunction(s) for diagram number 8619
    // (none)
    // Amplitude(s) for diagram number 8619
    FFV1_0( w_fp[3], w_fp[444], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];

    // *** DIAGRAM 8620 OF 15495 ***
    // Wavefunction(s) for diagram number 8620
    // (none)
    // Amplitude(s) for diagram number 8620
    FFV1_0( w_fp[523], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];

    // *** DIAGRAM 8621 OF 15495 ***
    // Wavefunction(s) for diagram number 8621
    // (none)
    // Amplitude(s) for diagram number 8621
    FFV1_0( w_fp[3], w_fp[244], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];

    // *** DIAGRAM 8622 OF 15495 ***
    // Wavefunction(s) for diagram number 8622
    // (none)
    // Amplitude(s) for diagram number 8622
    FFV1_0( w_fp[523], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8623 OF 15495 ***
    // Wavefunction(s) for diagram number 8623
    // (none)
    // Amplitude(s) for diagram number 8623
    FFV1_0( w_fp[184], w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];

    // *** DIAGRAM 8624 OF 15495 ***
    // Wavefunction(s) for diagram number 8624
    // (none)
    // Amplitude(s) for diagram number 8624
    FFV1_0( w_fp[184], w_fp[444], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8625 OF 15495 ***
    // Wavefunction(s) for diagram number 8625
    // (none)
    // Amplitude(s) for diagram number 8625
    FFV1_0( w_fp[3], w_fp[244], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[244], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[437] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[244], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];

    // *** DIAGRAM 8626 OF 15495 ***
    // Wavefunction(s) for diagram number 8626
    // (none)
    // Amplitude(s) for diagram number 8626
    FFV1_0( w_fp[184], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    FFV1_0( w_fp[184], w_fp[2], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    FFV1_0( w_fp[184], w_fp[2], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];

    // *** DIAGRAM 8627 OF 15495 ***
    // Wavefunction(s) for diagram number 8627
    // (none)
    // Amplitude(s) for diagram number 8627
    FFV1_0( w_fp[63], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[376], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[377], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8628 OF 15495 ***
    // Wavefunction(s) for diagram number 8628
    // (none)
    // Amplitude(s) for diagram number 8628
    FFV1_0( w_fp[3], w_fp[453], w_fp[303], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[172], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[76], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];

    // *** DIAGRAM 8629 OF 15495 ***
    // Wavefunction(s) for diagram number 8629
    // (none)
    // Amplitude(s) for diagram number 8629
    FFV1_0( w_fp[449], w_fp[484], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[513], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[446], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8630 OF 15495 ***
    // Wavefunction(s) for diagram number 8630
    // (none)
    // Amplitude(s) for diagram number 8630
    FFV1_0( w_fp[449], w_fp[2], w_fp[303], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[172], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[76], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 8631 OF 15495 ***
    // Wavefunction(s) for diagram number 8631
    // (none)
    // Amplitude(s) for diagram number 8631
    FFV1_0( w_fp[3], w_fp[484], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[513], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[175] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[446], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];

    // *** DIAGRAM 8632 OF 15495 ***
    // Wavefunction(s) for diagram number 8632
    // (none)
    // Amplitude(s) for diagram number 8632
    FFV1_0( w_fp[63], w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    FFV1_0( w_fp[376], w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    FFV1_0( w_fp[377], w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];

    // *** DIAGRAM 8633 OF 15495 ***
    // Wavefunction(s) for diagram number 8633
    // (none)
    // Amplitude(s) for diagram number 8633
    VVVV1_0( w_fp[532], w_fp[276], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[276], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[276], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[532], w_fp[356], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[356], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[356], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[532], w_fp[142], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[142], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[142], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8634 OF 15495 ***
    // Wavefunction(s) for diagram number 8634
    VVV1P0_1( w_fp[532], w_fp[276], COUPs[0], 1.0, depCoup, 0., 0., w_fp[559] );
    VVV1P0_1( w_fp[532], w_fp[356], COUPs[0], 1.0, depCoup, 0., 0., w_fp[560] );
    VVV1P0_1( w_fp[532], w_fp[142], COUPs[0], 1.0, depCoup, 0., 0., w_fp[536] );
    // Amplitude(s) for diagram number 8634
    VVV1_0( w_fp[9], w_fp[6], w_fp[559], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[560], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[536], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8635 OF 15495 ***
    // Wavefunction(s) for diagram number 8635
    // (none)
    // Amplitude(s) for diagram number 8635
    VVV1_0( w_fp[276], w_fp[6], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[356], w_fp[6], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[142], w_fp[6], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8636 OF 15495 ***
    // Wavefunction(s) for diagram number 8636
    // (none)
    // Amplitude(s) for diagram number 8636
    VVV1_0( w_fp[276], w_fp[9], w_fp[474], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[356], w_fp[9], w_fp[474], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[142], w_fp[9], w_fp[474], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8637 OF 15495 ***
    // Wavefunction(s) for diagram number 8637
    // (none)
    // Amplitude(s) for diagram number 8637
    FFV1_0( w_fp[3], w_fp[197], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[197], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];

    // *** DIAGRAM 8638 OF 15495 ***
    // Wavefunction(s) for diagram number 8638
    // (none)
    // Amplitude(s) for diagram number 8638
    FFV1_0( w_fp[3], w_fp[479], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[479], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[479], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8639 OF 15495 ***
    // Wavefunction(s) for diagram number 8639
    // (none)
    // Amplitude(s) for diagram number 8639
    FFV1_0( w_fp[449], w_fp[197], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[197], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[197], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8640 OF 15495 ***
    // Wavefunction(s) for diagram number 8640
    // (none)
    // Amplitude(s) for diagram number 8640
    FFV1_0( w_fp[174], w_fp[2], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[174], w_fp[2], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];

    // *** DIAGRAM 8641 OF 15495 ***
    // Wavefunction(s) for diagram number 8641
    // (none)
    // Amplitude(s) for diagram number 8641
    FFV1_0( w_fp[174], w_fp[453], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[453], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[453], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8642 OF 15495 ***
    // Wavefunction(s) for diagram number 8642
    // (none)
    // Amplitude(s) for diagram number 8642
    FFV1_0( w_fp[525], w_fp[2], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[525], w_fp[2], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[525], w_fp[2], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8643 OF 15495 ***
    // Wavefunction(s) for diagram number 8643
    // (none)
    // Amplitude(s) for diagram number 8643
    FFV1_0( w_fp[386], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[387], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[388], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8644 OF 15495 ***
    // Wavefunction(s) for diagram number 8644
    // (none)
    // Amplitude(s) for diagram number 8644
    FFV1_0( w_fp[3], w_fp[453], w_fp[378], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[381], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[382], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];

    // *** DIAGRAM 8645 OF 15495 ***
    // Wavefunction(s) for diagram number 8645
    // (none)
    // Amplitude(s) for diagram number 8645
    FFV1_0( w_fp[449], w_fp[458], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[440], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[455], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8646 OF 15495 ***
    // Wavefunction(s) for diagram number 8646
    // (none)
    // Amplitude(s) for diagram number 8646
    FFV1_0( w_fp[449], w_fp[2], w_fp[378], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[381], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[382], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];

    // *** DIAGRAM 8647 OF 15495 ***
    // Wavefunction(s) for diagram number 8647
    // (none)
    // Amplitude(s) for diagram number 8647
    FFV1_0( w_fp[3], w_fp[458], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[440], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[199] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[455], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];

    // *** DIAGRAM 8648 OF 15495 ***
    // Wavefunction(s) for diagram number 8648
    // (none)
    // Amplitude(s) for diagram number 8648
    FFV1_0( w_fp[386], w_fp[2], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    FFV1_0( w_fp[387], w_fp[2], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    FFV1_0( w_fp[388], w_fp[2], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];

    // *** DIAGRAM 8649 OF 15495 ***
    // Wavefunction(s) for diagram number 8649
    // (none)
    // Amplitude(s) for diagram number 8649
    VVVV1_0( w_fp[532], w_fp[93], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[93], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[93], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[532], w_fp[380], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[380], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[380], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[532], w_fp[379], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[379], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[379], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8650 OF 15495 ***
    // Wavefunction(s) for diagram number 8650
    VVV1P0_1( w_fp[532], w_fp[93], COUPs[0], 1.0, depCoup, 0., 0., w_fp[525] );
    VVV1P0_1( w_fp[532], w_fp[380], COUPs[0], 1.0, depCoup, 0., 0., w_fp[536] );
    VVV1P0_1( w_fp[532], w_fp[379], COUPs[0], 1.0, depCoup, 0., 0., w_fp[560] );
    // Amplitude(s) for diagram number 8650
    VVV1_0( w_fp[9], w_fp[5], w_fp[525], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[5], w_fp[536], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[5], w_fp[560], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8651 OF 15495 ***
    // Wavefunction(s) for diagram number 8651
    // (none)
    // Amplitude(s) for diagram number 8651
    VVV1_0( w_fp[93], w_fp[5], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[380], w_fp[5], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[379], w_fp[5], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8652 OF 15495 ***
    // Wavefunction(s) for diagram number 8652
    // (none)
    // Amplitude(s) for diagram number 8652
    VVV1_0( w_fp[93], w_fp[9], w_fp[505], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[380], w_fp[9], w_fp[505], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[379], w_fp[9], w_fp[505], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8653 OF 15495 ***
    // Wavefunction(s) for diagram number 8653
    // (none)
    // Amplitude(s) for diagram number 8653
    FFV1_0( w_fp[3], w_fp[169], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];

    // *** DIAGRAM 8654 OF 15495 ***
    // Wavefunction(s) for diagram number 8654
    // (none)
    // Amplitude(s) for diagram number 8654
    FFV1_0( w_fp[3], w_fp[588], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[588], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[588], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8655 OF 15495 ***
    // Wavefunction(s) for diagram number 8655
    // (none)
    // Amplitude(s) for diagram number 8655
    FFV1_0( w_fp[449], w_fp[169], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[169], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[169], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8656 OF 15495 ***
    // Wavefunction(s) for diagram number 8656
    // (none)
    // Amplitude(s) for diagram number 8656
    FFV1_0( w_fp[168], w_fp[2], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[168], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];

    // *** DIAGRAM 8657 OF 15495 ***
    // Wavefunction(s) for diagram number 8657
    // (none)
    // Amplitude(s) for diagram number 8657
    FFV1_0( w_fp[168], w_fp[453], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[453], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[453], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8658 OF 15495 ***
    // Wavefunction(s) for diagram number 8658
    // (none)
    // Amplitude(s) for diagram number 8658
    FFV1_0( w_fp[45], w_fp[2], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[45], w_fp[2], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[45], w_fp[2], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8659 OF 15495 ***
    // Wavefunction(s) for diagram number 8659
    // (none)
    // Amplitude(s) for diagram number 8659
    FFV1_0( w_fp[408], w_fp[453], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[409], w_fp[453], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[410], w_fp[453], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8660 OF 15495 ***
    // Wavefunction(s) for diagram number 8660
    // (none)
    // Amplitude(s) for diagram number 8660
    FFV1_0( w_fp[3], w_fp[453], w_fp[400], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[403], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[404], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];

    // *** DIAGRAM 8661 OF 15495 ***
    // Wavefunction(s) for diagram number 8661
    // (none)
    // Amplitude(s) for diagram number 8661
    FFV1_0( w_fp[449], w_fp[464], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[478], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[482], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8662 OF 15495 ***
    // Wavefunction(s) for diagram number 8662
    // (none)
    // Amplitude(s) for diagram number 8662
    FFV1_0( w_fp[449], w_fp[2], w_fp[400], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[403], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[404], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];

    // *** DIAGRAM 8663 OF 15495 ***
    // Wavefunction(s) for diagram number 8663
    // (none)
    // Amplitude(s) for diagram number 8663
    FFV1_0( w_fp[3], w_fp[464], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[478], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[205] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[482], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];

    // *** DIAGRAM 8664 OF 15495 ***
    // Wavefunction(s) for diagram number 8664
    // (none)
    // Amplitude(s) for diagram number 8664
    FFV1_0( w_fp[408], w_fp[2], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    FFV1_0( w_fp[409], w_fp[2], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    FFV1_0( w_fp[410], w_fp[2], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];

    // *** DIAGRAM 8665 OF 15495 ***
    // Wavefunction(s) for diagram number 8665
    // (none)
    // Amplitude(s) for diagram number 8665
    VVVV1_0( w_fp[532], w_fp[127], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[127], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[127], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[532], w_fp[402], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[402], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[402], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[532], w_fp[401], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[532], w_fp[401], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[532], w_fp[401], w_fp[9], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8666 OF 15495 ***
    // Wavefunction(s) for diagram number 8666
    VVV1P0_1( w_fp[532], w_fp[127], COUPs[0], 1.0, depCoup, 0., 0., w_fp[45] );
    VVV1P0_1( w_fp[532], w_fp[402], COUPs[0], 1.0, depCoup, 0., 0., w_fp[560] );
    VVV1P0_1( w_fp[532], w_fp[401], COUPs[0], 1.0, depCoup, 0., 0., w_fp[536] );
    // Amplitude(s) for diagram number 8666
    VVV1_0( w_fp[9], w_fp[4], w_fp[45], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[4], w_fp[560], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[4], w_fp[536], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8667 OF 15495 ***
    // Wavefunction(s) for diagram number 8667
    // (none)
    // Amplitude(s) for diagram number 8667
    VVV1_0( w_fp[127], w_fp[4], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[402], w_fp[4], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[401], w_fp[4], w_fp[541], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8668 OF 15495 ***
    // Wavefunction(s) for diagram number 8668
    // (none)
    // Amplitude(s) for diagram number 8668
    VVV1_0( w_fp[127], w_fp[9], w_fp[537], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[402], w_fp[9], w_fp[537], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[401], w_fp[9], w_fp[537], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8669 OF 15495 ***
    // Wavefunction(s) for diagram number 8669
    // (none)
    // Amplitude(s) for diagram number 8669
    FFV1_0( w_fp[3], w_fp[156], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[156], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];

    // *** DIAGRAM 8670 OF 15495 ***
    // Wavefunction(s) for diagram number 8670
    // (none)
    // Amplitude(s) for diagram number 8670
    FFV1_0( w_fp[3], w_fp[586], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[586], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[586], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8671 OF 15495 ***
    // Wavefunction(s) for diagram number 8671
    // (none)
    // Amplitude(s) for diagram number 8671
    FFV1_0( w_fp[449], w_fp[156], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[156], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[156], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8672 OF 15495 ***
    // Wavefunction(s) for diagram number 8672
    // (none)
    // Amplitude(s) for diagram number 8672
    FFV1_0( w_fp[196], w_fp[2], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[196], w_fp[2], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];

    // *** DIAGRAM 8673 OF 15495 ***
    // Wavefunction(s) for diagram number 8673
    // (none)
    // Amplitude(s) for diagram number 8673
    FFV1_0( w_fp[196], w_fp[453], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[453], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[453], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8674 OF 15495 ***
    // Wavefunction(s) for diagram number 8674
    // (none)
    // Amplitude(s) for diagram number 8674
    FFV1_0( w_fp[442], w_fp[2], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8675 OF 15495 ***
    // Wavefunction(s) for diagram number 8675
    // (none)
    // Amplitude(s) for diagram number 8675
    FFV1_0( w_fp[3], w_fp[195], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[67], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[325] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[15], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];

    // *** DIAGRAM 8676 OF 15495 ***
    // Wavefunction(s) for diagram number 8676
    // (none)
    // Amplitude(s) for diagram number 8676
    FFV1_0( w_fp[50], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    FFV1_0( w_fp[49], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    FFV1_0( w_fp[48], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];

    // *** DIAGRAM 8677 OF 15495 ***
    // Wavefunction(s) for diagram number 8677
    // (none)
    // Amplitude(s) for diagram number 8677
    FFV1_0( w_fp[3], w_fp[453], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[433], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[453], w_fp[432], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];

    // *** DIAGRAM 8678 OF 15495 ***
    // Wavefunction(s) for diagram number 8678
    // (none)
    // Amplitude(s) for diagram number 8678
    FFV1_0( w_fp[50], w_fp[453], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[49], w_fp[453], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[48], w_fp[453], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8679 OF 15495 ***
    // Wavefunction(s) for diagram number 8679
    // (none)
    // Amplitude(s) for diagram number 8679
    FFV1_0( w_fp[449], w_fp[2], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[433], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    FFV1_0( w_fp[449], w_fp[2], w_fp[432], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 8680 OF 15495 ***
    // Wavefunction(s) for diagram number 8680
    // (none)
    // Amplitude(s) for diagram number 8680
    FFV1_0( w_fp[449], w_fp[195], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[67], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[449], w_fp[15], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8681 OF 15495 ***
    // Wavefunction(s) for diagram number 8681
    FFV1_1( w_fp[259], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[449] );
    FFV1_1( w_fp[449], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[453] );
    // Amplitude(s) for diagram number 8681
    FFV1_0( w_fp[216], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] -= amp_sv[0];

    // *** DIAGRAM 8682 OF 15495 ***
    // Wavefunction(s) for diagram number 8682
    FFV1_1( w_fp[449], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[585] );
    // Amplitude(s) for diagram number 8682
    FFV1_0( w_fp[216], w_fp[585], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 8683 OF 15495 ***
    // Wavefunction(s) for diagram number 8683
    FFV1_1( w_fp[449], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[442] );
    // Amplitude(s) for diagram number 8683
    FFV1_0( w_fp[198], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= amp_sv[0];

    // *** DIAGRAM 8684 OF 15495 ***
    // Wavefunction(s) for diagram number 8684
    // (none)
    // Amplitude(s) for diagram number 8684
    FFV1_0( w_fp[198], w_fp[585], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[141] -= amp_sv[0];

    // *** DIAGRAM 8685 OF 15495 ***
    // Wavefunction(s) for diagram number 8685
    // (none)
    // Amplitude(s) for diagram number 8685
    FFV1_0( w_fp[199], w_fp[442], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= amp_sv[0];

    // *** DIAGRAM 8686 OF 15495 ***
    // Wavefunction(s) for diagram number 8686
    // (none)
    // Amplitude(s) for diagram number 8686
    FFV1_0( w_fp[199], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] -= amp_sv[0];

    // *** DIAGRAM 8687 OF 15495 ***
    // Wavefunction(s) for diagram number 8687
    FFV1_2( w_fp[196], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[536] );
    FFV1_2( w_fp[536], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[560] );
    // Amplitude(s) for diagram number 8687
    FFV1_0( w_fp[560], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= amp_sv[0];

    // *** DIAGRAM 8688 OF 15495 ***
    // Wavefunction(s) for diagram number 8688
    FFV1_2( w_fp[536], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[45] );
    // Amplitude(s) for diagram number 8688
    FFV1_0( w_fp[45], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];

    // *** DIAGRAM 8689 OF 15495 ***
    // Wavefunction(s) for diagram number 8689
    FFV1_2( w_fp[536], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[586] );
    // Amplitude(s) for diagram number 8689
    FFV1_0( w_fp[586], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] -= amp_sv[0];

    // *** DIAGRAM 8690 OF 15495 ***
    // Wavefunction(s) for diagram number 8690
    // (none)
    // Amplitude(s) for diagram number 8690
    FFV1_0( w_fp[45], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] -= amp_sv[0];

    // *** DIAGRAM 8691 OF 15495 ***
    // Wavefunction(s) for diagram number 8691
    // (none)
    // Amplitude(s) for diagram number 8691
    FFV1_0( w_fp[586], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[238] -= amp_sv[0];

    // *** DIAGRAM 8692 OF 15495 ***
    // Wavefunction(s) for diagram number 8692
    // (none)
    // Amplitude(s) for diagram number 8692
    FFV1_0( w_fp[560], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 8693 OF 15495 ***
    // Wavefunction(s) for diagram number 8693
    FFV1_1( w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[537] );
    // Amplitude(s) for diagram number 8693
    FFV1_0( w_fp[198], w_fp[537], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[173] -= amp_sv[0];

    // *** DIAGRAM 8694 OF 15495 ***
    // Wavefunction(s) for diagram number 8694
    FFV1_2( w_fp[198], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[541] );
    // Amplitude(s) for diagram number 8694
    FFV1_0( w_fp[541], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] -= amp_sv[0];

    // *** DIAGRAM 8695 OF 15495 ***
    // Wavefunction(s) for diagram number 8695
    // (none)
    // Amplitude(s) for diagram number 8695
    FFV1_0( w_fp[199], w_fp[537], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] -= amp_sv[0];

    // *** DIAGRAM 8696 OF 15495 ***
    // Wavefunction(s) for diagram number 8696
    FFV1_2( w_fp[199], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[532] );
    // Amplitude(s) for diagram number 8696
    FFV1_0( w_fp[532], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] -= amp_sv[0];

    // *** DIAGRAM 8697 OF 15495 ***
    // Wavefunction(s) for diagram number 8697
    FFV1_1( w_fp[439], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[525] );
    // Amplitude(s) for diagram number 8697
    FFV1_0( w_fp[216], w_fp[525], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[197] -= amp_sv[0];

    // *** DIAGRAM 8698 OF 15495 ***
    // Wavefunction(s) for diagram number 8698
    FFV1_2( w_fp[216], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[588] );
    // Amplitude(s) for diagram number 8698
    FFV1_0( w_fp[588], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= amp_sv[0];

    // *** DIAGRAM 8699 OF 15495 ***
    // Wavefunction(s) for diagram number 8699
    // (none)
    // Amplitude(s) for diagram number 8699
    FFV1_0( w_fp[199], w_fp[525], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] -= amp_sv[0];

    // *** DIAGRAM 8700 OF 15495 ***
    // Wavefunction(s) for diagram number 8700
    // (none)
    // Amplitude(s) for diagram number 8700
    FFV1_0( w_fp[532], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[205] -= amp_sv[0];

    // *** DIAGRAM 8701 OF 15495 ***
    // Wavefunction(s) for diagram number 8701
    FFV1_1( w_fp[441], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[505] );
    // Amplitude(s) for diagram number 8701
    FFV1_0( w_fp[216], w_fp[505], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] -= amp_sv[0];

    // *** DIAGRAM 8702 OF 15495 ***
    // Wavefunction(s) for diagram number 8702
    // (none)
    // Amplitude(s) for diagram number 8702
    FFV1_0( w_fp[588], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 8703 OF 15495 ***
    // Wavefunction(s) for diagram number 8703
    // (none)
    // Amplitude(s) for diagram number 8703
    FFV1_0( w_fp[198], w_fp[505], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] -= amp_sv[0];

    // *** DIAGRAM 8704 OF 15495 ***
    // Wavefunction(s) for diagram number 8704
    // (none)
    // Amplitude(s) for diagram number 8704
    FFV1_0( w_fp[541], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] -= amp_sv[0];

    // *** DIAGRAM 8705 OF 15495 ***
    // Wavefunction(s) for diagram number 8705
    // (none)
    // Amplitude(s) for diagram number 8705
    FFV1_0( w_fp[243], w_fp[449], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8706 OF 15495 ***
    // Wavefunction(s) for diagram number 8706
    // (none)
    // Amplitude(s) for diagram number 8706
    FFV1_0( w_fp[199], w_fp[449], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8707 OF 15495 ***
    // Wavefunction(s) for diagram number 8707
    // (none)
    // Amplitude(s) for diagram number 8707
    FFV1_0( w_fp[196], w_fp[449], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 8708 OF 15495 ***
    // Wavefunction(s) for diagram number 8708
    // (none)
    // Amplitude(s) for diagram number 8708
    FFV1_0( w_fp[536], w_fp[465], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8709 OF 15495 ***
    // Wavefunction(s) for diagram number 8709
    // (none)
    // Amplitude(s) for diagram number 8709
    FFV1_0( w_fp[536], w_fp[441], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8710 OF 15495 ***
    // Wavefunction(s) for diagram number 8710
    // (none)
    // Amplitude(s) for diagram number 8710
    FFV1_0( w_fp[536], w_fp[259], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 8711 OF 15495 ***
    // Wavefunction(s) for diagram number 8711
    VVV1P0_1( w_fp[0], w_fp[113], COUPs[0], 1.0, depCoup, 0., 0., w_fp[559] );
    // Amplitude(s) for diagram number 8711
    VVV1_0( w_fp[559], w_fp[549], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8712 OF 15495 ***
    // Wavefunction(s) for diagram number 8712
    // (none)
    // Amplitude(s) for diagram number 8712
    FFV1_0( w_fp[196], w_fp[441], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 8713 OF 15495 ***
    // Wavefunction(s) for diagram number 8713
    // (none)
    // Amplitude(s) for diagram number 8713
    FFV1_0( w_fp[199], w_fp[259], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];

    // *** DIAGRAM 8714 OF 15495 ***
    // Wavefunction(s) for diagram number 8714
    // (none)
    // Amplitude(s) for diagram number 8714
    VVV1_0( w_fp[0], w_fp[549], w_fp[116], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8715 OF 15495 ***
    // Wavefunction(s) for diagram number 8715
    // (none)
    // Amplitude(s) for diagram number 8715
    FFV1_0( w_fp[199], w_fp[465], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8716 OF 15495 ***
    // Wavefunction(s) for diagram number 8716
    // (none)
    // Amplitude(s) for diagram number 8716
    FFV1_0( w_fp[243], w_fp[441], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8717 OF 15495 ***
    // Wavefunction(s) for diagram number 8717
    VVVV1P0_1( w_fp[0], w_fp[113], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[479] );
    VVVV3P0_1( w_fp[0], w_fp[113], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[474] );
    VVVV4P0_1( w_fp[0], w_fp[113], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[444] );
    // Amplitude(s) for diagram number 8717
    FFV1_0( w_fp[196], w_fp[259], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8718 OF 15495 ***
    // Wavefunction(s) for diagram number 8718
    // (none)
    // Amplitude(s) for diagram number 8718
    FFV1_0( w_fp[123], w_fp[449], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8719 OF 15495 ***
    // Wavefunction(s) for diagram number 8719
    // (none)
    // Amplitude(s) for diagram number 8719
    FFV1_0( w_fp[198], w_fp[449], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8720 OF 15495 ***
    // Wavefunction(s) for diagram number 8720
    // (none)
    // Amplitude(s) for diagram number 8720
    FFV1_0( w_fp[196], w_fp[449], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];

    // *** DIAGRAM 8721 OF 15495 ***
    // Wavefunction(s) for diagram number 8721
    // (none)
    // Amplitude(s) for diagram number 8721
    FFV1_0( w_fp[536], w_fp[466], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8722 OF 15495 ***
    // Wavefunction(s) for diagram number 8722
    // (none)
    // Amplitude(s) for diagram number 8722
    FFV1_0( w_fp[536], w_fp[439], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8723 OF 15495 ***
    // Wavefunction(s) for diagram number 8723
    // (none)
    // Amplitude(s) for diagram number 8723
    FFV1_0( w_fp[536], w_fp[259], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 8724 OF 15495 ***
    // Wavefunction(s) for diagram number 8724
    VVV1P0_1( w_fp[0], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[522] );
    // Amplitude(s) for diagram number 8724
    VVV1_0( w_fp[522], w_fp[549], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8725 OF 15495 ***
    // Wavefunction(s) for diagram number 8725
    // (none)
    // Amplitude(s) for diagram number 8725
    FFV1_0( w_fp[196], w_fp[439], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];

    // *** DIAGRAM 8726 OF 15495 ***
    // Wavefunction(s) for diagram number 8726
    // (none)
    // Amplitude(s) for diagram number 8726
    FFV1_0( w_fp[198], w_fp[259], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];

    // *** DIAGRAM 8727 OF 15495 ***
    // Wavefunction(s) for diagram number 8727
    // (none)
    // Amplitude(s) for diagram number 8727
    VVV1_0( w_fp[0], w_fp[549], w_fp[125], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8728 OF 15495 ***
    // Wavefunction(s) for diagram number 8728
    // (none)
    // Amplitude(s) for diagram number 8728
    FFV1_0( w_fp[198], w_fp[466], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8729 OF 15495 ***
    // Wavefunction(s) for diagram number 8729
    // (none)
    // Amplitude(s) for diagram number 8729
    FFV1_0( w_fp[123], w_fp[439], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8730 OF 15495 ***
    // Wavefunction(s) for diagram number 8730
    VVVV1P0_1( w_fp[0], w_fp[100], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[523] );
    VVVV3P0_1( w_fp[0], w_fp[100], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[580] );
    VVVV4P0_1( w_fp[0], w_fp[100], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[528] );
    // Amplitude(s) for diagram number 8730
    FFV1_0( w_fp[196], w_fp[259], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8731 OF 15495 ***
    // Wavefunction(s) for diagram number 8731
    // (none)
    // Amplitude(s) for diagram number 8731
    FFV1_0( w_fp[216], w_fp[449], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8732 OF 15495 ***
    // Wavefunction(s) for diagram number 8732
    // (none)
    // Amplitude(s) for diagram number 8732
    FFV1_0( w_fp[129], w_fp[449], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8733 OF 15495 ***
    // Wavefunction(s) for diagram number 8733
    // (none)
    // Amplitude(s) for diagram number 8733
    FFV1_0( w_fp[196], w_fp[449], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 8734 OF 15495 ***
    // Wavefunction(s) for diagram number 8734
    // (none)
    // Amplitude(s) for diagram number 8734
    FFV1_0( w_fp[536], w_fp[436], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8735 OF 15495 ***
    // Wavefunction(s) for diagram number 8735
    // (none)
    // Amplitude(s) for diagram number 8735
    FFV1_0( w_fp[536], w_fp[467], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8736 OF 15495 ***
    // Wavefunction(s) for diagram number 8736
    // (none)
    // Amplitude(s) for diagram number 8736
    FFV1_0( w_fp[536], w_fp[259], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 8737 OF 15495 ***
    // Wavefunction(s) for diagram number 8737
    VVV1P0_1( w_fp[0], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[516] );
    // Amplitude(s) for diagram number 8737
    VVV1_0( w_fp[516], w_fp[549], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8738 OF 15495 ***
    // Wavefunction(s) for diagram number 8738
    // (none)
    // Amplitude(s) for diagram number 8738
    FFV1_0( w_fp[196], w_fp[436], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];

    // *** DIAGRAM 8739 OF 15495 ***
    // Wavefunction(s) for diagram number 8739
    // (none)
    // Amplitude(s) for diagram number 8739
    FFV1_0( w_fp[216], w_fp[259], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 8740 OF 15495 ***
    // Wavefunction(s) for diagram number 8740
    // (none)
    // Amplitude(s) for diagram number 8740
    VVV1_0( w_fp[0], w_fp[549], w_fp[131], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8741 OF 15495 ***
    // Wavefunction(s) for diagram number 8741
    // (none)
    // Amplitude(s) for diagram number 8741
    FFV1_0( w_fp[129], w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8742 OF 15495 ***
    // Wavefunction(s) for diagram number 8742
    // (none)
    // Amplitude(s) for diagram number 8742
    FFV1_0( w_fp[216], w_fp[467], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8743 OF 15495 ***
    // Wavefunction(s) for diagram number 8743
    VVVV1P0_1( w_fp[0], w_fp[5], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[549] );
    VVVV3P0_1( w_fp[0], w_fp[5], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[451] );
    VVVV4P0_1( w_fp[0], w_fp[5], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[571] );
    // Amplitude(s) for diagram number 8743
    FFV1_0( w_fp[196], w_fp[259], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8744 OF 15495 ***
    // Wavefunction(s) for diagram number 8744
    // (none)
    // Amplitude(s) for diagram number 8744
    FFV1_0( w_fp[196], w_fp[449], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[449], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[449], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 8745 OF 15495 ***
    // Wavefunction(s) for diagram number 8745
    // (none)
    // Amplitude(s) for diagram number 8745
    FFV1_0( w_fp[536], w_fp[259], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    FFV1_0( w_fp[536], w_fp[259], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    FFV1_0( w_fp[536], w_fp[259], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];

    // *** DIAGRAM 8746 OF 15495 ***
    // Wavefunction(s) for diagram number 8746
    VVV1P0_1( w_fp[0], w_fp[151], COUPs[0], 1.0, depCoup, 0., 0., w_fp[452] );
    VVV1P0_1( w_fp[0], w_fp[152], COUPs[0], 1.0, depCoup, 0., 0., w_fp[488] );
    VVV1P0_1( w_fp[0], w_fp[153], COUPs[0], 1.0, depCoup, 0., 0., w_fp[437] );
    // Amplitude(s) for diagram number 8746
    FFV1_0( w_fp[196], w_fp[259], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8747 OF 15495 ***
    // Wavefunction(s) for diagram number 8747
    // (none)
    // Amplitude(s) for diagram number 8747
    FFV1_0( w_fp[218], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] -= amp_sv[0];

    // *** DIAGRAM 8748 OF 15495 ***
    // Wavefunction(s) for diagram number 8748
    // (none)
    // Amplitude(s) for diagram number 8748
    FFV1_0( w_fp[218], w_fp[585], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[142] -= amp_sv[0];

    // *** DIAGRAM 8749 OF 15495 ***
    // Wavefunction(s) for diagram number 8749
    FFV1_1( w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[546] );
    // Amplitude(s) for diagram number 8749
    FFV1_0( w_fp[170], w_fp[546], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= amp_sv[0];

    // *** DIAGRAM 8750 OF 15495 ***
    // Wavefunction(s) for diagram number 8750
    // (none)
    // Amplitude(s) for diagram number 8750
    FFV1_0( w_fp[170], w_fp[585], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] -= amp_sv[0];

    // *** DIAGRAM 8751 OF 15495 ***
    // Wavefunction(s) for diagram number 8751
    // (none)
    // Amplitude(s) for diagram number 8751
    FFV1_0( w_fp[171], w_fp[546], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];

    // *** DIAGRAM 8752 OF 15495 ***
    // Wavefunction(s) for diagram number 8752
    // (none)
    // Amplitude(s) for diagram number 8752
    FFV1_0( w_fp[171], w_fp[453], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];

    // *** DIAGRAM 8753 OF 15495 ***
    // Wavefunction(s) for diagram number 8753
    FFV1_2( w_fp[168], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[527] );
    FFV1_2( w_fp[527], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[538] );
    // Amplitude(s) for diagram number 8753
    FFV1_0( w_fp[538], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= amp_sv[0];

    // *** DIAGRAM 8754 OF 15495 ***
    // Wavefunction(s) for diagram number 8754
    FFV1_2( w_fp[527], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[543] );
    // Amplitude(s) for diagram number 8754
    FFV1_0( w_fp[543], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= amp_sv[0];

    // *** DIAGRAM 8755 OF 15495 ***
    // Wavefunction(s) for diagram number 8755
    FFV1_2( w_fp[527], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[509] );
    // Amplitude(s) for diagram number 8755
    FFV1_0( w_fp[509], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= amp_sv[0];

    // *** DIAGRAM 8756 OF 15495 ***
    // Wavefunction(s) for diagram number 8756
    // (none)
    // Amplitude(s) for diagram number 8756
    FFV1_0( w_fp[543], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= amp_sv[0];

    // *** DIAGRAM 8757 OF 15495 ***
    // Wavefunction(s) for diagram number 8757
    // (none)
    // Amplitude(s) for diagram number 8757
    FFV1_0( w_fp[509], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 8758 OF 15495 ***
    // Wavefunction(s) for diagram number 8758
    // (none)
    // Amplitude(s) for diagram number 8758
    FFV1_0( w_fp[538], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 8759 OF 15495 ***
    // Wavefunction(s) for diagram number 8759
    FFV1_1( w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[520] );
    // Amplitude(s) for diagram number 8759
    FFV1_0( w_fp[170], w_fp[520], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[149] -= amp_sv[0];

    // *** DIAGRAM 8760 OF 15495 ***
    // Wavefunction(s) for diagram number 8760
    FFV1_2( w_fp[170], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[535] );
    // Amplitude(s) for diagram number 8760
    FFV1_0( w_fp[535], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] -= amp_sv[0];

    // *** DIAGRAM 8761 OF 15495 ***
    // Wavefunction(s) for diagram number 8761
    // (none)
    // Amplitude(s) for diagram number 8761
    FFV1_0( w_fp[171], w_fp[520], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] -= amp_sv[0];

    // *** DIAGRAM 8762 OF 15495 ***
    // Wavefunction(s) for diagram number 8762
    FFV1_2( w_fp[171], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[450] );
    // Amplitude(s) for diagram number 8762
    FFV1_0( w_fp[450], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] -= amp_sv[0];

    // *** DIAGRAM 8763 OF 15495 ***
    // Wavefunction(s) for diagram number 8763
    // (none)
    // Amplitude(s) for diagram number 8763
    FFV1_0( w_fp[218], w_fp[525], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[196] -= amp_sv[0];

    // *** DIAGRAM 8764 OF 15495 ***
    // Wavefunction(s) for diagram number 8764
    FFV1_2( w_fp[218], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[481] );
    // Amplitude(s) for diagram number 8764
    FFV1_0( w_fp[481], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];

    // *** DIAGRAM 8765 OF 15495 ***
    // Wavefunction(s) for diagram number 8765
    // (none)
    // Amplitude(s) for diagram number 8765
    FFV1_0( w_fp[171], w_fp[525], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] -= amp_sv[0];

    // *** DIAGRAM 8766 OF 15495 ***
    // Wavefunction(s) for diagram number 8766
    // (none)
    // Amplitude(s) for diagram number 8766
    FFV1_0( w_fp[450], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[199] -= amp_sv[0];

    // *** DIAGRAM 8767 OF 15495 ***
    // Wavefunction(s) for diagram number 8767
    // (none)
    // Amplitude(s) for diagram number 8767
    FFV1_0( w_fp[218], w_fp[505], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] -= amp_sv[0];

    // *** DIAGRAM 8768 OF 15495 ***
    // Wavefunction(s) for diagram number 8768
    // (none)
    // Amplitude(s) for diagram number 8768
    FFV1_0( w_fp[481], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= amp_sv[0];

    // *** DIAGRAM 8769 OF 15495 ***
    // Wavefunction(s) for diagram number 8769
    // (none)
    // Amplitude(s) for diagram number 8769
    FFV1_0( w_fp[170], w_fp[505], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= amp_sv[0];

    // *** DIAGRAM 8770 OF 15495 ***
    // Wavefunction(s) for diagram number 8770
    // (none)
    // Amplitude(s) for diagram number 8770
    FFV1_0( w_fp[535], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[223] -= amp_sv[0];

    // *** DIAGRAM 8771 OF 15495 ***
    // Wavefunction(s) for diagram number 8771
    // (none)
    // Amplitude(s) for diagram number 8771
    FFV1_0( w_fp[200], w_fp[449], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8772 OF 15495 ***
    // Wavefunction(s) for diagram number 8772
    // (none)
    // Amplitude(s) for diagram number 8772
    FFV1_0( w_fp[171], w_fp[449], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8773 OF 15495 ***
    // Wavefunction(s) for diagram number 8773
    // (none)
    // Amplitude(s) for diagram number 8773
    FFV1_0( w_fp[168], w_fp[449], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 8774 OF 15495 ***
    // Wavefunction(s) for diagram number 8774
    // (none)
    // Amplitude(s) for diagram number 8774
    FFV1_0( w_fp[527], w_fp[460], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8775 OF 15495 ***
    // Wavefunction(s) for diagram number 8775
    // (none)
    // Amplitude(s) for diagram number 8775
    FFV1_0( w_fp[527], w_fp[441], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8776 OF 15495 ***
    // Wavefunction(s) for diagram number 8776
    // (none)
    // Amplitude(s) for diagram number 8776
    FFV1_0( w_fp[527], w_fp[259], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 8777 OF 15495 ***
    // Wavefunction(s) for diagram number 8777
    VVV1P0_1( w_fp[0], w_fp[86], COUPs[0], 1.0, depCoup, 0., 0., w_fp[576] );
    // Amplitude(s) for diagram number 8777
    VVV1_0( w_fp[576], w_fp[468], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8778 OF 15495 ***
    // Wavefunction(s) for diagram number 8778
    // (none)
    // Amplitude(s) for diagram number 8778
    FFV1_0( w_fp[168], w_fp[441], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 8779 OF 15495 ***
    // Wavefunction(s) for diagram number 8779
    // (none)
    // Amplitude(s) for diagram number 8779
    FFV1_0( w_fp[171], w_fp[259], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];

    // *** DIAGRAM 8780 OF 15495 ***
    // Wavefunction(s) for diagram number 8780
    // (none)
    // Amplitude(s) for diagram number 8780
    VVV1_0( w_fp[0], w_fp[468], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8781 OF 15495 ***
    // Wavefunction(s) for diagram number 8781
    // (none)
    // Amplitude(s) for diagram number 8781
    FFV1_0( w_fp[171], w_fp[460], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8782 OF 15495 ***
    // Wavefunction(s) for diagram number 8782
    // (none)
    // Amplitude(s) for diagram number 8782
    FFV1_0( w_fp[200], w_fp[441], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8783 OF 15495 ***
    // Wavefunction(s) for diagram number 8783
    VVVV1P0_1( w_fp[0], w_fp[86], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[550] );
    VVVV3P0_1( w_fp[0], w_fp[86], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[480] );
    VVVV4P0_1( w_fp[0], w_fp[86], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[558] );
    // Amplitude(s) for diagram number 8783
    FFV1_0( w_fp[168], w_fp[259], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8784 OF 15495 ***
    // Wavefunction(s) for diagram number 8784
    // (none)
    // Amplitude(s) for diagram number 8784
    FFV1_0( w_fp[99], w_fp[449], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8785 OF 15495 ***
    // Wavefunction(s) for diagram number 8785
    // (none)
    // Amplitude(s) for diagram number 8785
    FFV1_0( w_fp[170], w_fp[449], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8786 OF 15495 ***
    // Wavefunction(s) for diagram number 8786
    // (none)
    // Amplitude(s) for diagram number 8786
    FFV1_0( w_fp[168], w_fp[449], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];

    // *** DIAGRAM 8787 OF 15495 ***
    // Wavefunction(s) for diagram number 8787
    // (none)
    // Amplitude(s) for diagram number 8787
    FFV1_0( w_fp[527], w_fp[463], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8788 OF 15495 ***
    // Wavefunction(s) for diagram number 8788
    // (none)
    // Amplitude(s) for diagram number 8788
    FFV1_0( w_fp[527], w_fp[439], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8789 OF 15495 ***
    // Wavefunction(s) for diagram number 8789
    // (none)
    // Amplitude(s) for diagram number 8789
    FFV1_0( w_fp[527], w_fp[259], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 8790 OF 15495 ***
    // Wavefunction(s) for diagram number 8790
    VVV1P0_1( w_fp[0], w_fp[102], COUPs[0], 1.0, depCoup, 0., 0., w_fp[518] );
    // Amplitude(s) for diagram number 8790
    VVV1_0( w_fp[518], w_fp[468], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8791 OF 15495 ***
    // Wavefunction(s) for diagram number 8791
    // (none)
    // Amplitude(s) for diagram number 8791
    FFV1_0( w_fp[168], w_fp[439], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];

    // *** DIAGRAM 8792 OF 15495 ***
    // Wavefunction(s) for diagram number 8792
    // (none)
    // Amplitude(s) for diagram number 8792
    FFV1_0( w_fp[170], w_fp[259], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];

    // *** DIAGRAM 8793 OF 15495 ***
    // Wavefunction(s) for diagram number 8793
    // (none)
    // Amplitude(s) for diagram number 8793
    VVV1_0( w_fp[0], w_fp[468], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8794 OF 15495 ***
    // Wavefunction(s) for diagram number 8794
    // (none)
    // Amplitude(s) for diagram number 8794
    FFV1_0( w_fp[170], w_fp[463], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8795 OF 15495 ***
    // Wavefunction(s) for diagram number 8795
    // (none)
    // Amplitude(s) for diagram number 8795
    FFV1_0( w_fp[99], w_fp[439], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8796 OF 15495 ***
    // Wavefunction(s) for diagram number 8796
    VVVV1P0_1( w_fp[0], w_fp[102], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[486] );
    VVVV3P0_1( w_fp[0], w_fp[102], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[547] );
    VVVV4P0_1( w_fp[0], w_fp[102], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[529] );
    // Amplitude(s) for diagram number 8796
    FFV1_0( w_fp[168], w_fp[259], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[529], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8797 OF 15495 ***
    // Wavefunction(s) for diagram number 8797
    // (none)
    // Amplitude(s) for diagram number 8797
    FFV1_0( w_fp[218], w_fp[449], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8798 OF 15495 ***
    // Wavefunction(s) for diagram number 8798
    // (none)
    // Amplitude(s) for diagram number 8798
    FFV1_0( w_fp[112], w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8799 OF 15495 ***
    // Wavefunction(s) for diagram number 8799
    // (none)
    // Amplitude(s) for diagram number 8799
    FFV1_0( w_fp[168], w_fp[449], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 8800 OF 15495 ***
    // Wavefunction(s) for diagram number 8800
    // (none)
    // Amplitude(s) for diagram number 8800
    FFV1_0( w_fp[527], w_fp[443], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 45 );
    storeWf( wfs, w_cx, nevt, 437 );
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 444 );
    storeWf( wfs, w_cx, nevt, 449 );
    storeWf( wfs, w_cx, nevt, 450 );
    storeWf( wfs, w_cx, nevt, 451 );
    storeWf( wfs, w_cx, nevt, 452 );
    storeWf( wfs, w_cx, nevt, 453 );
    storeWf( wfs, w_cx, nevt, 474 );
    storeWf( wfs, w_cx, nevt, 479 );
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 486 );
    storeWf( wfs, w_cx, nevt, 488 );
    storeWf( wfs, w_cx, nevt, 505 );
    storeWf( wfs, w_cx, nevt, 509 );
    storeWf( wfs, w_cx, nevt, 516 );
    storeWf( wfs, w_cx, nevt, 518 );
    storeWf( wfs, w_cx, nevt, 520 );
    storeWf( wfs, w_cx, nevt, 522 );
    storeWf( wfs, w_cx, nevt, 523 );
    storeWf( wfs, w_cx, nevt, 525 );
    storeWf( wfs, w_cx, nevt, 527 );
    storeWf( wfs, w_cx, nevt, 528 );
    storeWf( wfs, w_cx, nevt, 529 );
    storeWf( wfs, w_cx, nevt, 532 );
    storeWf( wfs, w_cx, nevt, 535 );
    storeWf( wfs, w_cx, nevt, 536 );
    storeWf( wfs, w_cx, nevt, 537 );
    storeWf( wfs, w_cx, nevt, 538 );
    storeWf( wfs, w_cx, nevt, 541 );
    storeWf( wfs, w_cx, nevt, 543 );
    storeWf( wfs, w_cx, nevt, 546 );
    storeWf( wfs, w_cx, nevt, 547 );
    storeWf( wfs, w_cx, nevt, 549 );
    storeWf( wfs, w_cx, nevt, 550 );
    storeWf( wfs, w_cx, nevt, 558 );
    storeWf( wfs, w_cx, nevt, 559 );
    storeWf( wfs, w_cx, nevt, 560 );
    storeWf( wfs, w_cx, nevt, 571 );
    storeWf( wfs, w_cx, nevt, 576 );
    storeWf( wfs, w_cx, nevt, 580 );
    storeWf( wfs, w_cx, nevt, 585 );
    storeWf( wfs, w_cx, nevt, 586 );
    storeWf( wfs, w_cx, nevt, 588 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
