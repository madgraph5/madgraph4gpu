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
  diagramgroup42( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 164 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 173 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 178 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 185 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 201 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 203 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 207 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 217 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 219 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 222 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 512 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 596 );
#endif
#endif

    // *** DIAGRAM 8201 OF 15495 ***
    // Wavefunction(s) for diagram number 8201
    // (none)
    // Amplitude(s) for diagram number 8201
    FFV1_0( w_fp[450], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8202 OF 15495 ***
    // Wavefunction(s) for diagram number 8202
    // (none)
    // Amplitude(s) for diagram number 8202
    FFV1_0( w_fp[535], w_fp[512], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[277] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];

    // *** DIAGRAM 8203 OF 15495 ***
    // Wavefunction(s) for diagram number 8203
    // (none)
    // Amplitude(s) for diagram number 8203
    FFV1_0( w_fp[3], w_fp[512], w_fp[557], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8204 OF 15495 ***
    // Wavefunction(s) for diagram number 8204
    // (none)
    // Amplitude(s) for diagram number 8204
    VVVV1_0( w_fp[505], w_fp[1], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    VVVV3_0( w_fp[505], w_fp[1], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    VVVV4_0( w_fp[505], w_fp[1], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];

    // *** DIAGRAM 8205 OF 15495 ***
    // Wavefunction(s) for diagram number 8205
    // (none)
    // Amplitude(s) for diagram number 8205
    VVV1_0( w_fp[164], w_fp[6], w_fp[520], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];

    // *** DIAGRAM 8206 OF 15495 ***
    // Wavefunction(s) for diagram number 8206
    // (none)
    // Amplitude(s) for diagram number 8206
    VVV1_0( w_fp[1], w_fp[164], w_fp[557], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];

    // *** DIAGRAM 8207 OF 15495 ***
    // Wavefunction(s) for diagram number 8207
    // (none)
    // Amplitude(s) for diagram number 8207
    FFV1_0( w_fp[3], w_fp[161], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8208 OF 15495 ***
    // Wavefunction(s) for diagram number 8208
    // (none)
    // Amplitude(s) for diagram number 8208
    FFV1_0( w_fp[535], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[319] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];

    // *** DIAGRAM 8209 OF 15495 ***
    // Wavefunction(s) for diagram number 8209
    // (none)
    // Amplitude(s) for diagram number 8209
    FFV1_0( w_fp[547], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[271] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];

    // *** DIAGRAM 8210 OF 15495 ***
    // Wavefunction(s) for diagram number 8210
    // (none)
    // Amplitude(s) for diagram number 8210
    FFV1_0( w_fp[3], w_fp[512], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8211 OF 15495 ***
    // Wavefunction(s) for diagram number 8211
    // (none)
    // Amplitude(s) for diagram number 8211
    VVVV1_0( w_fp[474], w_fp[1], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[1], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[1], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];

    // *** DIAGRAM 8212 OF 15495 ***
    // Wavefunction(s) for diagram number 8212
    // (none)
    // Amplitude(s) for diagram number 8212
    VVV1_0( w_fp[164], w_fp[5], w_fp[486], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];

    // *** DIAGRAM 8213 OF 15495 ***
    // Wavefunction(s) for diagram number 8213
    // (none)
    // Amplitude(s) for diagram number 8213
    VVV1_0( w_fp[1], w_fp[164], w_fp[435], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];

    // *** DIAGRAM 8214 OF 15495 ***
    // Wavefunction(s) for diagram number 8214
    // (none)
    // Amplitude(s) for diagram number 8214
    FFV1_0( w_fp[3], w_fp[158], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8215 OF 15495 ***
    // Wavefunction(s) for diagram number 8215
    // (none)
    // Amplitude(s) for diagram number 8215
    FFV1_0( w_fp[547], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[295] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];

    // *** DIAGRAM 8216 OF 15495 ***
    // Wavefunction(s) for diagram number 8216
    // (none)
    // Amplitude(s) for diagram number 8216
    VVV1_0( w_fp[543], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    VVV1_0( w_fp[538], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    VVV1_0( w_fp[527], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];

    // *** DIAGRAM 8217 OF 15495 ***
    // Wavefunction(s) for diagram number 8217
    // (none)
    // Amplitude(s) for diagram number 8217
    FFV1_0( w_fp[3], w_fp[161], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[161], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[161], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8218 OF 15495 ***
    // Wavefunction(s) for diagram number 8218
    // (none)
    // Amplitude(s) for diagram number 8218
    VVV1_0( w_fp[558], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    VVV1_0( w_fp[480], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    VVV1_0( w_fp[550], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];

    // *** DIAGRAM 8219 OF 15495 ***
    // Wavefunction(s) for diagram number 8219
    // (none)
    // Amplitude(s) for diagram number 8219
    FFV1_0( w_fp[3], w_fp[158], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[158], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[158], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8220 OF 15495 ***
    // Wavefunction(s) for diagram number 8220
    // (none)
    // Amplitude(s) for diagram number 8220
    FFV1_0( w_fp[3], w_fp[512], w_fp[495], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[512], w_fp[596], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8221 OF 15495 ***
    // Wavefunction(s) for diagram number 8221
    // (none)
    // Amplitude(s) for diagram number 8221
    VVV1_0( w_fp[495], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    VVV1_0( w_fp[438], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[263] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    VVV1_0( w_fp[596], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];

    // *** DIAGRAM 8222 OF 15495 ***
    // Wavefunction(s) for diagram number 8222
    // (none)
    // Amplitude(s) for diagram number 8222
    VVV1_0( w_fp[585], w_fp[173], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8223 OF 15495 ***
    // Wavefunction(s) for diagram number 8223
    // (none)
    // Amplitude(s) for diagram number 8223
    FFV1_0( w_fp[168], w_fp[161], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];

    // *** DIAGRAM 8224 OF 15495 ***
    // Wavefunction(s) for diagram number 8224
    // (none)
    // Amplitude(s) for diagram number 8224
    FFV1_0( w_fp[170], w_fp[156], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];

    // *** DIAGRAM 8225 OF 15495 ***
    // Wavefunction(s) for diagram number 8225
    // (none)
    // Amplitude(s) for diagram number 8225
    FFV1_0( w_fp[256], w_fp[586], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8226 OF 15495 ***
    // Wavefunction(s) for diagram number 8226
    // (none)
    // Amplitude(s) for diagram number 8226
    FFV1_0( w_fp[170], w_fp[586], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8227 OF 15495 ***
    // Wavefunction(s) for diagram number 8227
    // (none)
    // Amplitude(s) for diagram number 8227
    FFV1_0( w_fp[45], w_fp[512], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8228 OF 15495 ***
    // Wavefunction(s) for diagram number 8228
    // (none)
    // Amplitude(s) for diagram number 8228
    FFV1_0( w_fp[45], w_fp[161], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8229 OF 15495 ***
    // Wavefunction(s) for diagram number 8229
    // (none)
    // Amplitude(s) for diagram number 8229
    FFV1_0( w_fp[168], w_fp[512], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[269] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];

    // *** DIAGRAM 8230 OF 15495 ***
    // Wavefunction(s) for diagram number 8230
    // (none)
    // Amplitude(s) for diagram number 8230
    FFV1_0( w_fp[256], w_fp[156], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];

    // *** DIAGRAM 8231 OF 15495 ***
    // Wavefunction(s) for diagram number 8231
    // (none)
    // Amplitude(s) for diagram number 8231
    VVV1_0( w_fp[474], w_fp[1], w_fp[173], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8232 OF 15495 ***
    // Wavefunction(s) for diagram number 8232
    // (none)
    // Amplitude(s) for diagram number 8232
    FFV1_0( w_fp[170], w_fp[512], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8233 OF 15495 ***
    // Wavefunction(s) for diagram number 8233
    // (none)
    // Amplitude(s) for diagram number 8233
    FFV1_0( w_fp[256], w_fp[161], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8234 OF 15495 ***
    // Wavefunction(s) for diagram number 8234
    // (none)
    // Amplitude(s) for diagram number 8234
    FFV1_0( w_fp[168], w_fp[156], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[156], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[156], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8235 OF 15495 ***
    // Wavefunction(s) for diagram number 8235
    // (none)
    // Amplitude(s) for diagram number 8235
    VVV1_0( w_fp[585], w_fp[178], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8236 OF 15495 ***
    // Wavefunction(s) for diagram number 8236
    // (none)
    // Amplitude(s) for diagram number 8236
    FFV1_0( w_fp[174], w_fp[158], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];

    // *** DIAGRAM 8237 OF 15495 ***
    // Wavefunction(s) for diagram number 8237
    // (none)
    // Amplitude(s) for diagram number 8237
    FFV1_0( w_fp[175], w_fp[156], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];

    // *** DIAGRAM 8238 OF 15495 ***
    // Wavefunction(s) for diagram number 8238
    // (none)
    // Amplitude(s) for diagram number 8238
    FFV1_0( w_fp[255], w_fp[586], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8239 OF 15495 ***
    // Wavefunction(s) for diagram number 8239
    // (none)
    // Amplitude(s) for diagram number 8239
    FFV1_0( w_fp[175], w_fp[586], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8240 OF 15495 ***
    // Wavefunction(s) for diagram number 8240
    // (none)
    // Amplitude(s) for diagram number 8240
    FFV1_0( w_fp[525], w_fp[512], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8241 OF 15495 ***
    // Wavefunction(s) for diagram number 8241
    // (none)
    // Amplitude(s) for diagram number 8241
    FFV1_0( w_fp[525], w_fp[158], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8242 OF 15495 ***
    // Wavefunction(s) for diagram number 8242
    // (none)
    // Amplitude(s) for diagram number 8242
    FFV1_0( w_fp[174], w_fp[512], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];

    // *** DIAGRAM 8243 OF 15495 ***
    // Wavefunction(s) for diagram number 8243
    // (none)
    // Amplitude(s) for diagram number 8243
    FFV1_0( w_fp[255], w_fp[156], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];

    // *** DIAGRAM 8244 OF 15495 ***
    // Wavefunction(s) for diagram number 8244
    // (none)
    // Amplitude(s) for diagram number 8244
    VVV1_0( w_fp[505], w_fp[1], w_fp[178], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8245 OF 15495 ***
    // Wavefunction(s) for diagram number 8245
    // (none)
    // Amplitude(s) for diagram number 8245
    FFV1_0( w_fp[175], w_fp[512], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8246 OF 15495 ***
    // Wavefunction(s) for diagram number 8246
    // (none)
    // Amplitude(s) for diagram number 8246
    FFV1_0( w_fp[255], w_fp[158], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8247 OF 15495 ***
    // Wavefunction(s) for diagram number 8247
    // (none)
    // Amplitude(s) for diagram number 8247
    FFV1_0( w_fp[174], w_fp[156], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[156], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[156], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8248 OF 15495 ***
    // Wavefunction(s) for diagram number 8248
    // (none)
    // Amplitude(s) for diagram number 8248
    VVV1_0( w_fp[585], w_fp[164], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];

    // *** DIAGRAM 8249 OF 15495 ***
    // Wavefunction(s) for diagram number 8249
    // (none)
    // Amplitude(s) for diagram number 8249
    FFV1_0( w_fp[3], w_fp[185], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8250 OF 15495 ***
    // Wavefunction(s) for diagram number 8250
    // (none)
    // Amplitude(s) for diagram number 8250
    FFV1_0( w_fp[184], w_fp[156], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8251 OF 15495 ***
    // Wavefunction(s) for diagram number 8251
    // (none)
    // Amplitude(s) for diagram number 8251
    FFV1_0( w_fp[3], w_fp[586], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8252 OF 15495 ***
    // Wavefunction(s) for diagram number 8252
    // (none)
    // Amplitude(s) for diagram number 8252
    FFV1_0( w_fp[184], w_fp[586], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];

    // *** DIAGRAM 8253 OF 15495 ***
    // Wavefunction(s) for diagram number 8253
    // (none)
    // Amplitude(s) for diagram number 8253
    FFV1_0( w_fp[449], w_fp[512], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[272] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];

    // *** DIAGRAM 8254 OF 15495 ***
    // Wavefunction(s) for diagram number 8254
    // (none)
    // Amplitude(s) for diagram number 8254
    FFV1_0( w_fp[449], w_fp[156], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8255 OF 15495 ***
    // Wavefunction(s) for diagram number 8255
    // (none)
    // Amplitude(s) for diagram number 8255
    FFV1_0( w_fp[449], w_fp[185], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];

    // *** DIAGRAM 8256 OF 15495 ***
    // Wavefunction(s) for diagram number 8256
    // (none)
    // Amplitude(s) for diagram number 8256
    FFV1_0( w_fp[3], w_fp[512], w_fp[531], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8257 OF 15495 ***
    // Wavefunction(s) for diagram number 8257
    // (none)
    // Amplitude(s) for diagram number 8257
    VVV1_0( w_fp[531], w_fp[1], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];

    // *** DIAGRAM 8258 OF 15495 ***
    // Wavefunction(s) for diagram number 8258
    // (none)
    // Amplitude(s) for diagram number 8258
    FFV1_0( w_fp[184], w_fp[512], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];

    // *** DIAGRAM 8259 OF 15495 ***
    // Wavefunction(s) for diagram number 8259
    // (none)
    // Amplitude(s) for diagram number 8259
    VVV1_0( w_fp[532], w_fp[359], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8260 OF 15495 ***
    // Wavefunction(s) for diagram number 8260
    // (none)
    // Amplitude(s) for diagram number 8260
    FFV1_0( w_fp[3], w_fp[156], w_fp[539], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[156], w_fp[591], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[156], w_fp[592], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];

    // *** DIAGRAM 8261 OF 15495 ***
    // Wavefunction(s) for diagram number 8261
    // (none)
    // Amplitude(s) for diagram number 8261
    VVVV1_0( w_fp[585], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    VVVV3_0( w_fp[585], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    VVVV4_0( w_fp[585], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];

    // *** DIAGRAM 8262 OF 15495 ***
    // Wavefunction(s) for diagram number 8262
    // (none)
    // Amplitude(s) for diagram number 8262
    VVV1_0( w_fp[194], w_fp[6], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];

    // *** DIAGRAM 8263 OF 15495 ***
    // Wavefunction(s) for diagram number 8263
    // (none)
    // Amplitude(s) for diagram number 8263
    VVV1_0( w_fp[194], w_fp[4], w_fp[556], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];

    // *** DIAGRAM 8264 OF 15495 ***
    // Wavefunction(s) for diagram number 8264
    // (none)
    // Amplitude(s) for diagram number 8264
    FFV1_0( w_fp[488], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[421] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];

    // *** DIAGRAM 8265 OF 15495 ***
    // Wavefunction(s) for diagram number 8265
    // (none)
    // Amplitude(s) for diagram number 8265
    FFV1_0( w_fp[3], w_fp[190], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8266 OF 15495 ***
    // Wavefunction(s) for diagram number 8266
    // (none)
    // Amplitude(s) for diagram number 8266
    FFV1_0( w_fp[488], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];

    // *** DIAGRAM 8267 OF 15495 ***
    // Wavefunction(s) for diagram number 8267
    // (none)
    // Amplitude(s) for diagram number 8267
    FFV1_0( w_fp[3], w_fp[191], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8268 OF 15495 ***
    // Wavefunction(s) for diagram number 8268
    // (none)
    // Amplitude(s) for diagram number 8268
    FFV1_0( w_fp[516], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8269 OF 15495 ***
    // Wavefunction(s) for diagram number 8269
    // (none)
    // Amplitude(s) for diagram number 8269
    FFV1_0( w_fp[529], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8270 OF 15495 ***
    // Wavefunction(s) for diagram number 8270
    // (none)
    // Amplitude(s) for diagram number 8270
    FFV1_0( w_fp[528], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8271 OF 15495 ***
    // Wavefunction(s) for diagram number 8271
    // (none)
    // Amplitude(s) for diagram number 8271
    FFV1_0( w_fp[529], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8272 OF 15495 ***
    // Wavefunction(s) for diagram number 8272
    // (none)
    // Amplitude(s) for diagram number 8272
    FFV1_0( w_fp[528], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8273 OF 15495 ***
    // Wavefunction(s) for diagram number 8273
    // (none)
    // Amplitude(s) for diagram number 8273
    FFV1_0( w_fp[516], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8274 OF 15495 ***
    // Wavefunction(s) for diagram number 8274
    // (none)
    // Amplitude(s) for diagram number 8274
    FFV1_0( w_fp[523], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[397] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];

    // *** DIAGRAM 8275 OF 15495 ***
    // Wavefunction(s) for diagram number 8275
    // (none)
    // Amplitude(s) for diagram number 8275
    FFV1_0( w_fp[3], w_fp[477], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8276 OF 15495 ***
    // Wavefunction(s) for diagram number 8276
    // (none)
    // Amplitude(s) for diagram number 8276
    VVVV1_0( w_fp[537], w_fp[1], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    VVVV3_0( w_fp[537], w_fp[1], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    VVVV4_0( w_fp[537], w_fp[1], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];

    // *** DIAGRAM 8277 OF 15495 ***
    // Wavefunction(s) for diagram number 8277
    // (none)
    // Amplitude(s) for diagram number 8277
    VVV1_0( w_fp[194], w_fp[6], w_fp[522], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];

    // *** DIAGRAM 8278 OF 15495 ***
    // Wavefunction(s) for diagram number 8278
    // (none)
    // Amplitude(s) for diagram number 8278
    VVV1_0( w_fp[1], w_fp[194], w_fp[553], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];

    // *** DIAGRAM 8279 OF 15495 ***
    // Wavefunction(s) for diagram number 8279
    // (none)
    // Amplitude(s) for diagram number 8279
    FFV1_0( w_fp[3], w_fp[191], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8280 OF 15495 ***
    // Wavefunction(s) for diagram number 8280
    // (none)
    // Amplitude(s) for diagram number 8280
    FFV1_0( w_fp[523], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[439] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];

    // *** DIAGRAM 8281 OF 15495 ***
    // Wavefunction(s) for diagram number 8281
    // (none)
    // Amplitude(s) for diagram number 8281
    FFV1_0( w_fp[547], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[391] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];

    // *** DIAGRAM 8282 OF 15495 ***
    // Wavefunction(s) for diagram number 8282
    // (none)
    // Amplitude(s) for diagram number 8282
    FFV1_0( w_fp[3], w_fp[477], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8283 OF 15495 ***
    // Wavefunction(s) for diagram number 8283
    // (none)
    // Amplitude(s) for diagram number 8283
    VVVV1_0( w_fp[474], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];

    // *** DIAGRAM 8284 OF 15495 ***
    // Wavefunction(s) for diagram number 8284
    // (none)
    // Amplitude(s) for diagram number 8284
    VVV1_0( w_fp[194], w_fp[4], w_fp[486], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];

    // *** DIAGRAM 8285 OF 15495 ***
    // Wavefunction(s) for diagram number 8285
    // (none)
    // Amplitude(s) for diagram number 8285
    VVV1_0( w_fp[1], w_fp[194], w_fp[582], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];

    // *** DIAGRAM 8286 OF 15495 ***
    // Wavefunction(s) for diagram number 8286
    // (none)
    // Amplitude(s) for diagram number 8286
    FFV1_0( w_fp[3], w_fp[190], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8287 OF 15495 ***
    // Wavefunction(s) for diagram number 8287
    // (none)
    // Amplitude(s) for diagram number 8287
    FFV1_0( w_fp[547], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[415] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];

    // *** DIAGRAM 8288 OF 15495 ***
    // Wavefunction(s) for diagram number 8288
    // (none)
    // Amplitude(s) for diagram number 8288
    VVV1_0( w_fp[536], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    VVV1_0( w_fp[560], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    VVV1_0( w_fp[559], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];

    // *** DIAGRAM 8289 OF 15495 ***
    // Wavefunction(s) for diagram number 8289
    // (none)
    // Amplitude(s) for diagram number 8289
    FFV1_0( w_fp[3], w_fp[191], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[191], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[191], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8290 OF 15495 ***
    // Wavefunction(s) for diagram number 8290
    // (none)
    // Amplitude(s) for diagram number 8290
    VVV1_0( w_fp[558], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    VVV1_0( w_fp[480], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    VVV1_0( w_fp[550], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];

    // *** DIAGRAM 8291 OF 15495 ***
    // Wavefunction(s) for diagram number 8291
    // (none)
    // Amplitude(s) for diagram number 8291
    FFV1_0( w_fp[3], w_fp[190], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[190], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[190], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8292 OF 15495 ***
    // Wavefunction(s) for diagram number 8292
    // (none)
    // Amplitude(s) for diagram number 8292
    FFV1_0( w_fp[3], w_fp[477], w_fp[471], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8293 OF 15495 ***
    // Wavefunction(s) for diagram number 8293
    // (none)
    // Amplitude(s) for diagram number 8293
    VVV1_0( w_fp[471], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    VVV1_0( w_fp[562], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[383] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    VVV1_0( w_fp[524], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];

    // *** DIAGRAM 8294 OF 15495 ***
    // Wavefunction(s) for diagram number 8294
    // (none)
    // Amplitude(s) for diagram number 8294
    VVV1_0( w_fp[585], w_fp[201], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8295 OF 15495 ***
    // Wavefunction(s) for diagram number 8295
    // (none)
    // Amplitude(s) for diagram number 8295
    FFV1_0( w_fp[196], w_fp[191], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];

    // *** DIAGRAM 8296 OF 15495 ***
    // Wavefunction(s) for diagram number 8296
    // (none)
    // Amplitude(s) for diagram number 8296
    FFV1_0( w_fp[198], w_fp[169], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];

    // *** DIAGRAM 8297 OF 15495 ***
    // Wavefunction(s) for diagram number 8297
    // (none)
    // Amplitude(s) for diagram number 8297
    FFV1_0( w_fp[355], w_fp[588], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8298 OF 15495 ***
    // Wavefunction(s) for diagram number 8298
    // (none)
    // Amplitude(s) for diagram number 8298
    FFV1_0( w_fp[198], w_fp[588], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8299 OF 15495 ***
    // Wavefunction(s) for diagram number 8299
    // (none)
    // Amplitude(s) for diagram number 8299
    FFV1_0( w_fp[442], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8300 OF 15495 ***
    // Wavefunction(s) for diagram number 8300
    // (none)
    // Amplitude(s) for diagram number 8300
    FFV1_0( w_fp[442], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8301 OF 15495 ***
    // Wavefunction(s) for diagram number 8301
    // (none)
    // Amplitude(s) for diagram number 8301
    FFV1_0( w_fp[196], w_fp[477], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[389] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];

    // *** DIAGRAM 8302 OF 15495 ***
    // Wavefunction(s) for diagram number 8302
    // (none)
    // Amplitude(s) for diagram number 8302
    FFV1_0( w_fp[355], w_fp[169], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];

    // *** DIAGRAM 8303 OF 15495 ***
    // Wavefunction(s) for diagram number 8303
    // (none)
    // Amplitude(s) for diagram number 8303
    VVV1_0( w_fp[474], w_fp[1], w_fp[201], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8304 OF 15495 ***
    // Wavefunction(s) for diagram number 8304
    // (none)
    // Amplitude(s) for diagram number 8304
    FFV1_0( w_fp[198], w_fp[477], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8305 OF 15495 ***
    // Wavefunction(s) for diagram number 8305
    // (none)
    // Amplitude(s) for diagram number 8305
    FFV1_0( w_fp[355], w_fp[191], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8306 OF 15495 ***
    // Wavefunction(s) for diagram number 8306
    // (none)
    // Amplitude(s) for diagram number 8306
    FFV1_0( w_fp[196], w_fp[169], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[169], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[169], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8307 OF 15495 ***
    // Wavefunction(s) for diagram number 8307
    // (none)
    // Amplitude(s) for diagram number 8307
    VVV1_0( w_fp[585], w_fp[203], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8308 OF 15495 ***
    // Wavefunction(s) for diagram number 8308
    // (none)
    // Amplitude(s) for diagram number 8308
    FFV1_0( w_fp[174], w_fp[190], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];

    // *** DIAGRAM 8309 OF 15495 ***
    // Wavefunction(s) for diagram number 8309
    // (none)
    // Amplitude(s) for diagram number 8309
    FFV1_0( w_fp[202], w_fp[169], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];

    // *** DIAGRAM 8310 OF 15495 ***
    // Wavefunction(s) for diagram number 8310
    // (none)
    // Amplitude(s) for diagram number 8310
    FFV1_0( w_fp[255], w_fp[588], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8311 OF 15495 ***
    // Wavefunction(s) for diagram number 8311
    // (none)
    // Amplitude(s) for diagram number 8311
    FFV1_0( w_fp[202], w_fp[588], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8312 OF 15495 ***
    // Wavefunction(s) for diagram number 8312
    // (none)
    // Amplitude(s) for diagram number 8312
    FFV1_0( w_fp[525], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8313 OF 15495 ***
    // Wavefunction(s) for diagram number 8313
    // (none)
    // Amplitude(s) for diagram number 8313
    FFV1_0( w_fp[525], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8314 OF 15495 ***
    // Wavefunction(s) for diagram number 8314
    // (none)
    // Amplitude(s) for diagram number 8314
    FFV1_0( w_fp[174], w_fp[477], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];

    // *** DIAGRAM 8315 OF 15495 ***
    // Wavefunction(s) for diagram number 8315
    // (none)
    // Amplitude(s) for diagram number 8315
    FFV1_0( w_fp[255], w_fp[169], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];

    // *** DIAGRAM 8316 OF 15495 ***
    // Wavefunction(s) for diagram number 8316
    // (none)
    // Amplitude(s) for diagram number 8316
    VVV1_0( w_fp[537], w_fp[1], w_fp[203], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8317 OF 15495 ***
    // Wavefunction(s) for diagram number 8317
    // (none)
    // Amplitude(s) for diagram number 8317
    FFV1_0( w_fp[202], w_fp[477], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8318 OF 15495 ***
    // Wavefunction(s) for diagram number 8318
    // (none)
    // Amplitude(s) for diagram number 8318
    FFV1_0( w_fp[255], w_fp[190], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8319 OF 15495 ***
    // Wavefunction(s) for diagram number 8319
    // (none)
    // Amplitude(s) for diagram number 8319
    FFV1_0( w_fp[174], w_fp[169], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8320 OF 15495 ***
    // Wavefunction(s) for diagram number 8320
    // (none)
    // Amplitude(s) for diagram number 8320
    VVV1_0( w_fp[585], w_fp[194], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];

    // *** DIAGRAM 8321 OF 15495 ***
    // Wavefunction(s) for diagram number 8321
    // (none)
    // Amplitude(s) for diagram number 8321
    FFV1_0( w_fp[3], w_fp[207], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8322 OF 15495 ***
    // Wavefunction(s) for diagram number 8322
    // (none)
    // Amplitude(s) for diagram number 8322
    FFV1_0( w_fp[206], w_fp[169], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8323 OF 15495 ***
    // Wavefunction(s) for diagram number 8323
    // (none)
    // Amplitude(s) for diagram number 8323
    FFV1_0( w_fp[3], w_fp[588], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8324 OF 15495 ***
    // Wavefunction(s) for diagram number 8324
    // (none)
    // Amplitude(s) for diagram number 8324
    FFV1_0( w_fp[206], w_fp[588], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];

    // *** DIAGRAM 8325 OF 15495 ***
    // Wavefunction(s) for diagram number 8325
    // (none)
    // Amplitude(s) for diagram number 8325
    FFV1_0( w_fp[449], w_fp[477], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];

    // *** DIAGRAM 8326 OF 15495 ***
    // Wavefunction(s) for diagram number 8326
    // (none)
    // Amplitude(s) for diagram number 8326
    FFV1_0( w_fp[449], w_fp[169], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8327 OF 15495 ***
    // Wavefunction(s) for diagram number 8327
    // (none)
    // Amplitude(s) for diagram number 8327
    FFV1_0( w_fp[449], w_fp[207], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];

    // *** DIAGRAM 8328 OF 15495 ***
    // Wavefunction(s) for diagram number 8328
    // (none)
    // Amplitude(s) for diagram number 8328
    FFV1_0( w_fp[3], w_fp[477], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8329 OF 15495 ***
    // Wavefunction(s) for diagram number 8329
    // (none)
    // Amplitude(s) for diagram number 8329
    VVV1_0( w_fp[526], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];

    // *** DIAGRAM 8330 OF 15495 ***
    // Wavefunction(s) for diagram number 8330
    // (none)
    // Amplitude(s) for diagram number 8330
    FFV1_0( w_fp[206], w_fp[477], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];

    // *** DIAGRAM 8331 OF 15495 ***
    // Wavefunction(s) for diagram number 8331
    // (none)
    // Amplitude(s) for diagram number 8331
    VVV1_0( w_fp[532], w_fp[362], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8332 OF 15495 ***
    // Wavefunction(s) for diagram number 8332
    // (none)
    // Amplitude(s) for diagram number 8332
    FFV1_0( w_fp[3], w_fp[169], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[169], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];

    // *** DIAGRAM 8333 OF 15495 ***
    // Wavefunction(s) for diagram number 8333
    // (none)
    // Amplitude(s) for diagram number 8333
    VVVV1_0( w_fp[585], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    VVVV3_0( w_fp[585], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    VVVV4_0( w_fp[585], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];

    // *** DIAGRAM 8334 OF 15495 ***
    // Wavefunction(s) for diagram number 8334
    // (none)
    // Amplitude(s) for diagram number 8334
    VVV1_0( w_fp[214], w_fp[5], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];

    // *** DIAGRAM 8335 OF 15495 ***
    // Wavefunction(s) for diagram number 8335
    // (none)
    // Amplitude(s) for diagram number 8335
    VVV1_0( w_fp[214], w_fp[4], w_fp[576], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];

    // *** DIAGRAM 8336 OF 15495 ***
    // Wavefunction(s) for diagram number 8336
    // (none)
    // Amplitude(s) for diagram number 8336
    FFV1_0( w_fp[488], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[541] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];

    // *** DIAGRAM 8337 OF 15495 ***
    // Wavefunction(s) for diagram number 8337
    // (none)
    // Amplitude(s) for diagram number 8337
    FFV1_0( w_fp[3], w_fp[211], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8338 OF 15495 ***
    // Wavefunction(s) for diagram number 8338
    // (none)
    // Amplitude(s) for diagram number 8338
    FFV1_0( w_fp[488], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];

    // *** DIAGRAM 8339 OF 15495 ***
    // Wavefunction(s) for diagram number 8339
    // (none)
    // Amplitude(s) for diagram number 8339
    FFV1_0( w_fp[3], w_fp[212], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8340 OF 15495 ***
    // Wavefunction(s) for diagram number 8340
    // (none)
    // Amplitude(s) for diagram number 8340
    FFV1_0( w_fp[516], w_fp[476], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8341 OF 15495 ***
    // Wavefunction(s) for diagram number 8341
    // (none)
    // Amplitude(s) for diagram number 8341
    FFV1_0( w_fp[450], w_fp[476], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8342 OF 15495 ***
    // Wavefunction(s) for diagram number 8342
    // (none)
    // Amplitude(s) for diagram number 8342
    FFV1_0( w_fp[528], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8343 OF 15495 ***
    // Wavefunction(s) for diagram number 8343
    // (none)
    // Amplitude(s) for diagram number 8343
    FFV1_0( w_fp[450], w_fp[211], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8344 OF 15495 ***
    // Wavefunction(s) for diagram number 8344
    // (none)
    // Amplitude(s) for diagram number 8344
    FFV1_0( w_fp[528], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8345 OF 15495 ***
    // Wavefunction(s) for diagram number 8345
    // (none)
    // Amplitude(s) for diagram number 8345
    FFV1_0( w_fp[516], w_fp[212], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8346 OF 15495 ***
    // Wavefunction(s) for diagram number 8346
    // (none)
    // Amplitude(s) for diagram number 8346
    FFV1_0( w_fp[523], w_fp[476], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[517] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];

    // *** DIAGRAM 8347 OF 15495 ***
    // Wavefunction(s) for diagram number 8347
    // (none)
    // Amplitude(s) for diagram number 8347
    FFV1_0( w_fp[3], w_fp[476], w_fp[555], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8348 OF 15495 ***
    // Wavefunction(s) for diagram number 8348
    // (none)
    // Amplitude(s) for diagram number 8348
    VVVV1_0( w_fp[537], w_fp[1], w_fp[214], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    VVVV3_0( w_fp[537], w_fp[1], w_fp[214], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    VVVV4_0( w_fp[537], w_fp[1], w_fp[214], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];

    // *** DIAGRAM 8349 OF 15495 ***
    // Wavefunction(s) for diagram number 8349
    // (none)
    // Amplitude(s) for diagram number 8349
    VVV1_0( w_fp[214], w_fp[5], w_fp[522], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];

    // *** DIAGRAM 8350 OF 15495 ***
    // Wavefunction(s) for diagram number 8350
    // (none)
    // Amplitude(s) for diagram number 8350
    VVV1_0( w_fp[1], w_fp[214], w_fp[555], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];

    // *** DIAGRAM 8351 OF 15495 ***
    // Wavefunction(s) for diagram number 8351
    // (none)
    // Amplitude(s) for diagram number 8351
    FFV1_0( w_fp[3], w_fp[212], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8352 OF 15495 ***
    // Wavefunction(s) for diagram number 8352
    // (none)
    // Amplitude(s) for diagram number 8352
    FFV1_0( w_fp[523], w_fp[212], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[559] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];

    // *** DIAGRAM 8353 OF 15495 ***
    // Wavefunction(s) for diagram number 8353
    // (none)
    // Amplitude(s) for diagram number 8353
    FFV1_0( w_fp[535], w_fp[476], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[511] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];

    // *** DIAGRAM 8354 OF 15495 ***
    // Wavefunction(s) for diagram number 8354
    // (none)
    // Amplitude(s) for diagram number 8354
    FFV1_0( w_fp[3], w_fp[476], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8355 OF 15495 ***
    // Wavefunction(s) for diagram number 8355
    // (none)
    // Amplitude(s) for diagram number 8355
    VVVV1_0( w_fp[505], w_fp[1], w_fp[214], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    VVVV3_0( w_fp[505], w_fp[1], w_fp[214], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    VVVV4_0( w_fp[505], w_fp[1], w_fp[214], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 8356 OF 15495 ***
    // Wavefunction(s) for diagram number 8356
    // (none)
    // Amplitude(s) for diagram number 8356
    VVV1_0( w_fp[214], w_fp[4], w_fp[520], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];

    // *** DIAGRAM 8357 OF 15495 ***
    // Wavefunction(s) for diagram number 8357
    // (none)
    // Amplitude(s) for diagram number 8357
    VVV1_0( w_fp[1], w_fp[214], w_fp[530], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 8358 OF 15495 ***
    // Wavefunction(s) for diagram number 8358
    // (none)
    // Amplitude(s) for diagram number 8358
    FFV1_0( w_fp[3], w_fp[211], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8359 OF 15495 ***
    // Wavefunction(s) for diagram number 8359
    // (none)
    // Amplitude(s) for diagram number 8359
    FFV1_0( w_fp[535], w_fp[211], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[535] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];

    // *** DIAGRAM 8360 OF 15495 ***
    // Wavefunction(s) for diagram number 8360
    // (none)
    // Amplitude(s) for diagram number 8360
    VVV1_0( w_fp[536], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    VVV1_0( w_fp[560], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    VVV1_0( w_fp[559], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];

    // *** DIAGRAM 8361 OF 15495 ***
    // Wavefunction(s) for diagram number 8361
    // (none)
    // Amplitude(s) for diagram number 8361
    FFV1_0( w_fp[3], w_fp[212], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[212], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[212], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8362 OF 15495 ***
    // Wavefunction(s) for diagram number 8362
    // (none)
    // Amplitude(s) for diagram number 8362
    VVV1_0( w_fp[543], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    VVV1_0( w_fp[538], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    VVV1_0( w_fp[527], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];

    // *** DIAGRAM 8363 OF 15495 ***
    // Wavefunction(s) for diagram number 8363
    // (none)
    // Amplitude(s) for diagram number 8363
    FFV1_0( w_fp[3], w_fp[211], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[211], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[211], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8364 OF 15495 ***
    // Wavefunction(s) for diagram number 8364
    // (none)
    // Amplitude(s) for diagram number 8364
    FFV1_0( w_fp[3], w_fp[476], w_fp[595], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[476], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[476], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8365 OF 15495 ***
    // Wavefunction(s) for diagram number 8365
    // (none)
    // Amplitude(s) for diagram number 8365
    VVV1_0( w_fp[595], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    VVV1_0( w_fp[578], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[503] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    VVV1_0( w_fp[577], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];

    // *** DIAGRAM 8366 OF 15495 ***
    // Wavefunction(s) for diagram number 8366
    // (none)
    // Amplitude(s) for diagram number 8366
    VVV1_0( w_fp[585], w_fp[217], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8367 OF 15495 ***
    // Wavefunction(s) for diagram number 8367
    // (none)
    // Amplitude(s) for diagram number 8367
    FFV1_0( w_fp[196], w_fp[212], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];

    // *** DIAGRAM 8368 OF 15495 ***
    // Wavefunction(s) for diagram number 8368
    // (none)
    // Amplitude(s) for diagram number 8368
    FFV1_0( w_fp[216], w_fp[197], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];

    // *** DIAGRAM 8369 OF 15495 ***
    // Wavefunction(s) for diagram number 8369
    // (none)
    // Amplitude(s) for diagram number 8369
    FFV1_0( w_fp[355], w_fp[479], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8370 OF 15495 ***
    // Wavefunction(s) for diagram number 8370
    // (none)
    // Amplitude(s) for diagram number 8370
    FFV1_0( w_fp[216], w_fp[479], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8371 OF 15495 ***
    // Wavefunction(s) for diagram number 8371
    // (none)
    // Amplitude(s) for diagram number 8371
    FFV1_0( w_fp[442], w_fp[476], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8372 OF 15495 ***
    // Wavefunction(s) for diagram number 8372
    // (none)
    // Amplitude(s) for diagram number 8372
    FFV1_0( w_fp[442], w_fp[212], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8373 OF 15495 ***
    // Wavefunction(s) for diagram number 8373
    // (none)
    // Amplitude(s) for diagram number 8373
    FFV1_0( w_fp[196], w_fp[476], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[509] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];

    // *** DIAGRAM 8374 OF 15495 ***
    // Wavefunction(s) for diagram number 8374
    // (none)
    // Amplitude(s) for diagram number 8374
    FFV1_0( w_fp[355], w_fp[197], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];

    // *** DIAGRAM 8375 OF 15495 ***
    // Wavefunction(s) for diagram number 8375
    // (none)
    // Amplitude(s) for diagram number 8375
    VVV1_0( w_fp[505], w_fp[1], w_fp[217], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8376 OF 15495 ***
    // Wavefunction(s) for diagram number 8376
    // (none)
    // Amplitude(s) for diagram number 8376
    FFV1_0( w_fp[216], w_fp[476], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8377 OF 15495 ***
    // Wavefunction(s) for diagram number 8377
    // (none)
    // Amplitude(s) for diagram number 8377
    FFV1_0( w_fp[355], w_fp[212], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8378 OF 15495 ***
    // Wavefunction(s) for diagram number 8378
    // (none)
    // Amplitude(s) for diagram number 8378
    FFV1_0( w_fp[196], w_fp[197], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[197], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[197], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8379 OF 15495 ***
    // Wavefunction(s) for diagram number 8379
    // (none)
    // Amplitude(s) for diagram number 8379
    VVV1_0( w_fp[585], w_fp[219], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8380 OF 15495 ***
    // Wavefunction(s) for diagram number 8380
    // (none)
    // Amplitude(s) for diagram number 8380
    FFV1_0( w_fp[168], w_fp[211], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];

    // *** DIAGRAM 8381 OF 15495 ***
    // Wavefunction(s) for diagram number 8381
    // (none)
    // Amplitude(s) for diagram number 8381
    FFV1_0( w_fp[218], w_fp[197], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];

    // *** DIAGRAM 8382 OF 15495 ***
    // Wavefunction(s) for diagram number 8382
    // (none)
    // Amplitude(s) for diagram number 8382
    FFV1_0( w_fp[256], w_fp[479], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8383 OF 15495 ***
    // Wavefunction(s) for diagram number 8383
    // (none)
    // Amplitude(s) for diagram number 8383
    FFV1_0( w_fp[218], w_fp[479], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8384 OF 15495 ***
    // Wavefunction(s) for diagram number 8384
    // (none)
    // Amplitude(s) for diagram number 8384
    FFV1_0( w_fp[45], w_fp[476], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8385 OF 15495 ***
    // Wavefunction(s) for diagram number 8385
    // (none)
    // Amplitude(s) for diagram number 8385
    FFV1_0( w_fp[45], w_fp[211], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8386 OF 15495 ***
    // Wavefunction(s) for diagram number 8386
    // (none)
    // Amplitude(s) for diagram number 8386
    FFV1_0( w_fp[168], w_fp[476], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];

    // *** DIAGRAM 8387 OF 15495 ***
    // Wavefunction(s) for diagram number 8387
    // (none)
    // Amplitude(s) for diagram number 8387
    FFV1_0( w_fp[256], w_fp[197], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];

    // *** DIAGRAM 8388 OF 15495 ***
    // Wavefunction(s) for diagram number 8388
    // (none)
    // Amplitude(s) for diagram number 8388
    VVV1_0( w_fp[537], w_fp[1], w_fp[219], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8389 OF 15495 ***
    // Wavefunction(s) for diagram number 8389
    // (none)
    // Amplitude(s) for diagram number 8389
    FFV1_0( w_fp[218], w_fp[476], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8390 OF 15495 ***
    // Wavefunction(s) for diagram number 8390
    // (none)
    // Amplitude(s) for diagram number 8390
    FFV1_0( w_fp[256], w_fp[211], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8391 OF 15495 ***
    // Wavefunction(s) for diagram number 8391
    // (none)
    // Amplitude(s) for diagram number 8391
    FFV1_0( w_fp[168], w_fp[197], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[197], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[197], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8392 OF 15495 ***
    // Wavefunction(s) for diagram number 8392
    // (none)
    // Amplitude(s) for diagram number 8392
    VVV1_0( w_fp[585], w_fp[214], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];

    // *** DIAGRAM 8393 OF 15495 ***
    // Wavefunction(s) for diagram number 8393
    // (none)
    // Amplitude(s) for diagram number 8393
    FFV1_0( w_fp[3], w_fp[222], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8394 OF 15495 ***
    // Wavefunction(s) for diagram number 8394
    // (none)
    // Amplitude(s) for diagram number 8394
    FFV1_0( w_fp[221], w_fp[197], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8395 OF 15495 ***
    // Wavefunction(s) for diagram number 8395
    // (none)
    // Amplitude(s) for diagram number 8395
    FFV1_0( w_fp[3], w_fp[479], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8396 OF 15495 ***
    // Wavefunction(s) for diagram number 8396
    // (none)
    // Amplitude(s) for diagram number 8396
    FFV1_0( w_fp[221], w_fp[479], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];

    // *** DIAGRAM 8397 OF 15495 ***
    // Wavefunction(s) for diagram number 8397
    // (none)
    // Amplitude(s) for diagram number 8397
    FFV1_0( w_fp[449], w_fp[476], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];

    // *** DIAGRAM 8398 OF 15495 ***
    // Wavefunction(s) for diagram number 8398
    // (none)
    // Amplitude(s) for diagram number 8398
    FFV1_0( w_fp[449], w_fp[197], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 8399 OF 15495 ***
    // Wavefunction(s) for diagram number 8399
    // (none)
    // Amplitude(s) for diagram number 8399
    FFV1_0( w_fp[449], w_fp[222], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 8400 OF 15495 ***
    // Wavefunction(s) for diagram number 8400
    // (none)
    // Amplitude(s) for diagram number 8400
    FFV1_0( w_fp[3], w_fp[476], w_fp[594], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];

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
