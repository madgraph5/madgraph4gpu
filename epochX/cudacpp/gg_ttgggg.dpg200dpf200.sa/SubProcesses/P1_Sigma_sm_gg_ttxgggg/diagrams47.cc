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
  diagramgroup47( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 105 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 121 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 185 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 189 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 207 );
    retrieveWf( wfs, w_cx, nevt, 209 );
    retrieveWf( wfs, w_cx, nevt, 210 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 222 );
    retrieveWf( wfs, w_cx, nevt, 223 );
    retrieveWf( wfs, w_cx, nevt, 224 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 593 );
#endif
#endif

    // *** DIAGRAM 9201 OF 15495 ***
    // Wavefunction(s) for diagram number 9201
    // (none)
    // Amplitude(s) for diagram number 9201
    FFV1_0( w_fp[37], w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    FFV1_0( w_fp[40], w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    FFV1_0( w_fp[12], w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];

    // *** DIAGRAM 9202 OF 15495 ***
    // Wavefunction(s) for diagram number 9202
    VVVV1P0_1( w_fp[0], w_fp[145], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[645] );
    VVVV3P0_1( w_fp[0], w_fp[145], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[646] );
    VVVV4P0_1( w_fp[0], w_fp[145], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[647] );
    VVVV1P0_1( w_fp[0], w_fp[146], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[648] );
    VVVV3P0_1( w_fp[0], w_fp[146], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[649] );
    VVVV4P0_1( w_fp[0], w_fp[146], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[650] );
    VVVV1P0_1( w_fp[0], w_fp[147], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[651] );
    VVVV3P0_1( w_fp[0], w_fp[147], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[652] );
    VVVV4P0_1( w_fp[0], w_fp[147], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[653] );
    // Amplitude(s) for diagram number 9202
    FFV1_0( w_fp[3], w_fp[259], w_fp[645], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[646], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[647], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[648], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[649], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[650], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[651], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[652], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[653], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 9203 OF 15495 ***
    // Wavefunction(s) for diagram number 9203
    // (none)
    // Amplitude(s) for diagram number 9203
    FFV1_0( w_fp[26], w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    FFV1_0( w_fp[160], w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    FFV1_0( w_fp[105], w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];

    // *** DIAGRAM 9204 OF 15495 ***
    // Wavefunction(s) for diagram number 9204
    // (none)
    // Amplitude(s) for diagram number 9204
    FFV1_0( w_fp[3], w_fp[449], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[449], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[449], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9205 OF 15495 ***
    // Wavefunction(s) for diagram number 9205
    // (none)
    // Amplitude(s) for diagram number 9205
    VVV1_0( w_fp[452], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    VVV1_0( w_fp[488], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    VVV1_0( w_fp[437], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 9206 OF 15495 ***
    // Wavefunction(s) for diagram number 9206
    // (none)
    // Amplitude(s) for diagram number 9206
    FFV1_0( w_fp[3], w_fp[443], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[443], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[443], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9207 OF 15495 ***
    // Wavefunction(s) for diagram number 9207
    // (none)
    // Amplitude(s) for diagram number 9207
    VVV1_0( w_fp[0], w_fp[456], w_fp[58], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[456], w_fp[57], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[456], w_fp[38], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 9208 OF 15495 ***
    // Wavefunction(s) for diagram number 9208
    // (none)
    // Amplitude(s) for diagram number 9208
    FFV1_0( w_fp[26], w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    FFV1_0( w_fp[160], w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    FFV1_0( w_fp[105], w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];

    // *** DIAGRAM 9209 OF 15495 ***
    // Wavefunction(s) for diagram number 9209
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[151], COUPs[2], 1.0, depCoup, 0., 0., w_fp[449] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[151], COUPs[2], 1.0, depCoup, 0., 0., w_fp[654] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[151], COUPs[2], 1.0, depCoup, 0., 0., w_fp[655] );
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[152], COUPs[2], 1.0, depCoup, 0., 0., w_fp[656] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[152], COUPs[2], 1.0, depCoup, 0., 0., w_fp[657] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[152], COUPs[2], 1.0, depCoup, 0., 0., w_fp[658] );
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[153], COUPs[2], 1.0, depCoup, 0., 0., w_fp[659] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[153], COUPs[2], 1.0, depCoup, 0., 0., w_fp[660] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[153], COUPs[2], 1.0, depCoup, 0., 0., w_fp[661] );
    // Amplitude(s) for diagram number 9209
    FFV1_0( w_fp[3], w_fp[259], w_fp[449], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[654], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[655], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[656], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[657], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[658], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[659], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[660], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[661], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 9210 OF 15495 ***
    // Wavefunction(s) for diagram number 9210
    FFV1_2( w_fp[157], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[662] );
    FFV1_2( w_fp[662], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[663] );
    // Amplitude(s) for diagram number 9210
    FFV1_0( w_fp[663], w_fp[158], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= amp_sv[0];

    // *** DIAGRAM 9211 OF 15495 ***
    // Wavefunction(s) for diagram number 9211
    FFV1_2( w_fp[662], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[664] );
    // Amplitude(s) for diagram number 9211
    FFV1_0( w_fp[664], w_fp[158], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];

    // *** DIAGRAM 9212 OF 15495 ***
    // Wavefunction(s) for diagram number 9212
    FFV1_2( w_fp[662], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[665] );
    // Amplitude(s) for diagram number 9212
    FFV1_0( w_fp[665], w_fp[161], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= amp_sv[0];

    // *** DIAGRAM 9213 OF 15495 ***
    // Wavefunction(s) for diagram number 9213
    // (none)
    // Amplitude(s) for diagram number 9213
    FFV1_0( w_fp[664], w_fp[161], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= amp_sv[0];

    // *** DIAGRAM 9214 OF 15495 ***
    // Wavefunction(s) for diagram number 9214
    // (none)
    // Amplitude(s) for diagram number 9214
    FFV1_0( w_fp[665], w_fp[163], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] -= amp_sv[0];

    // *** DIAGRAM 9215 OF 15495 ***
    // Wavefunction(s) for diagram number 9215
    // (none)
    // Amplitude(s) for diagram number 9215
    FFV1_0( w_fp[663], w_fp[163], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= amp_sv[0];

    // *** DIAGRAM 9216 OF 15495 ***
    // Wavefunction(s) for diagram number 9216
    FFV1_1( w_fp[156], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[666] );
    FFV1_1( w_fp[666], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[667] );
    // Amplitude(s) for diagram number 9216
    FFV1_0( w_fp[29], w_fp[667], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] -= amp_sv[0];

    // *** DIAGRAM 9217 OF 15495 ***
    // Wavefunction(s) for diagram number 9217
    FFV1_1( w_fp[666], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[668] );
    // Amplitude(s) for diagram number 9217
    FFV1_0( w_fp[29], w_fp[668], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[263] -= amp_sv[0];

    // *** DIAGRAM 9218 OF 15495 ***
    // Wavefunction(s) for diagram number 9218
    FFV1_1( w_fp[666], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[669] );
    // Amplitude(s) for diagram number 9218
    FFV1_0( w_fp[39], w_fp[669], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];

    // *** DIAGRAM 9219 OF 15495 ***
    // Wavefunction(s) for diagram number 9219
    // (none)
    // Amplitude(s) for diagram number 9219
    FFV1_0( w_fp[39], w_fp[668], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] -= amp_sv[0];

    // *** DIAGRAM 9220 OF 15495 ***
    // Wavefunction(s) for diagram number 9220
    // (none)
    // Amplitude(s) for diagram number 9220
    FFV1_0( w_fp[77], w_fp[669], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];

    // *** DIAGRAM 9221 OF 15495 ***
    // Wavefunction(s) for diagram number 9221
    // (none)
    // Amplitude(s) for diagram number 9221
    FFV1_0( w_fp[77], w_fp[667], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] -= amp_sv[0];

    // *** DIAGRAM 9222 OF 15495 ***
    // Wavefunction(s) for diagram number 9222
    FFV1_2( w_fp[29], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[670] );
    // Amplitude(s) for diagram number 9222
    FFV1_0( w_fp[670], w_fp[161], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] -= amp_sv[0];

    // *** DIAGRAM 9223 OF 15495 ***
    // Wavefunction(s) for diagram number 9223
    FFV1_1( w_fp[161], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[671] );
    // Amplitude(s) for diagram number 9223
    FFV1_0( w_fp[29], w_fp[671], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[317] -= amp_sv[0];

    // *** DIAGRAM 9224 OF 15495 ***
    // Wavefunction(s) for diagram number 9224
    // (none)
    // Amplitude(s) for diagram number 9224
    FFV1_0( w_fp[670], w_fp[163], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[355] -= amp_sv[0];

    // *** DIAGRAM 9225 OF 15495 ***
    // Wavefunction(s) for diagram number 9225
    FFV1_1( w_fp[163], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[672] );
    // Amplitude(s) for diagram number 9225
    FFV1_0( w_fp[29], w_fp[672], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[341] -= amp_sv[0];

    // *** DIAGRAM 9226 OF 15495 ***
    // Wavefunction(s) for diagram number 9226
    FFV1_2( w_fp[39], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[673] );
    // Amplitude(s) for diagram number 9226
    FFV1_0( w_fp[673], w_fp[158], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] -= amp_sv[0];

    // *** DIAGRAM 9227 OF 15495 ***
    // Wavefunction(s) for diagram number 9227
    FFV1_1( w_fp[158], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[674] );
    // Amplitude(s) for diagram number 9227
    FFV1_0( w_fp[39], w_fp[674], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[293] -= amp_sv[0];

    // *** DIAGRAM 9228 OF 15495 ***
    // Wavefunction(s) for diagram number 9228
    // (none)
    // Amplitude(s) for diagram number 9228
    FFV1_0( w_fp[673], w_fp[163], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[349] -= amp_sv[0];

    // *** DIAGRAM 9229 OF 15495 ***
    // Wavefunction(s) for diagram number 9229
    // (none)
    // Amplitude(s) for diagram number 9229
    FFV1_0( w_fp[39], w_fp[672], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] -= amp_sv[0];

    // *** DIAGRAM 9230 OF 15495 ***
    // Wavefunction(s) for diagram number 9230
    FFV1_2( w_fp[77], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[675] );
    // Amplitude(s) for diagram number 9230
    FFV1_0( w_fp[675], w_fp[158], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] -= amp_sv[0];

    // *** DIAGRAM 9231 OF 15495 ***
    // Wavefunction(s) for diagram number 9231
    // (none)
    // Amplitude(s) for diagram number 9231
    FFV1_0( w_fp[77], w_fp[674], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] -= amp_sv[0];

    // *** DIAGRAM 9232 OF 15495 ***
    // Wavefunction(s) for diagram number 9232
    // (none)
    // Amplitude(s) for diagram number 9232
    FFV1_0( w_fp[675], w_fp[161], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[325] -= amp_sv[0];

    // *** DIAGRAM 9233 OF 15495 ***
    // Wavefunction(s) for diagram number 9233
    // (none)
    // Amplitude(s) for diagram number 9233
    FFV1_0( w_fp[77], w_fp[671], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] -= amp_sv[0];

    // *** DIAGRAM 9234 OF 15495 ***
    // Wavefunction(s) for diagram number 9234
    // (none)
    // Amplitude(s) for diagram number 9234
    FFV1_0( w_fp[662], w_fp[185], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9235 OF 15495 ***
    // Wavefunction(s) for diagram number 9235
    // (none)
    // Amplitude(s) for diagram number 9235
    FFV1_0( w_fp[662], w_fp[163], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9236 OF 15495 ***
    // Wavefunction(s) for diagram number 9236
    // (none)
    // Amplitude(s) for diagram number 9236
    FFV1_0( w_fp[662], w_fp[156], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 9237 OF 15495 ***
    // Wavefunction(s) for diagram number 9237
    // (none)
    // Amplitude(s) for diagram number 9237
    FFV1_0( w_fp[121], w_fp[666], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9238 OF 15495 ***
    // Wavefunction(s) for diagram number 9238
    // (none)
    // Amplitude(s) for diagram number 9238
    FFV1_0( w_fp[77], w_fp[666], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9239 OF 15495 ***
    // Wavefunction(s) for diagram number 9239
    // (none)
    // Amplitude(s) for diagram number 9239
    FFV1_0( w_fp[157], w_fp[666], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];

    // *** DIAGRAM 9240 OF 15495 ***
    // Wavefunction(s) for diagram number 9240
    // (none)
    // Amplitude(s) for diagram number 9240
    VVV1_0( w_fp[559], w_fp[497], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9241 OF 15495 ***
    // Wavefunction(s) for diagram number 9241
    // (none)
    // Amplitude(s) for diagram number 9241
    FFV1_0( w_fp[77], w_fp[156], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];

    // *** DIAGRAM 9242 OF 15495 ***
    // Wavefunction(s) for diagram number 9242
    // (none)
    // Amplitude(s) for diagram number 9242
    FFV1_0( w_fp[157], w_fp[163], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 9243 OF 15495 ***
    // Wavefunction(s) for diagram number 9243
    // (none)
    // Amplitude(s) for diagram number 9243
    VVV1_0( w_fp[0], w_fp[497], w_fp[116], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9244 OF 15495 ***
    // Wavefunction(s) for diagram number 9244
    // (none)
    // Amplitude(s) for diagram number 9244
    FFV1_0( w_fp[121], w_fp[163], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9245 OF 15495 ***
    // Wavefunction(s) for diagram number 9245
    // (none)
    // Amplitude(s) for diagram number 9245
    FFV1_0( w_fp[77], w_fp[185], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9246 OF 15495 ***
    // Wavefunction(s) for diagram number 9246
    // (none)
    // Amplitude(s) for diagram number 9246
    FFV1_0( w_fp[157], w_fp[156], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9247 OF 15495 ***
    // Wavefunction(s) for diagram number 9247
    // (none)
    // Amplitude(s) for diagram number 9247
    FFV1_0( w_fp[662], w_fp[187], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9248 OF 15495 ***
    // Wavefunction(s) for diagram number 9248
    // (none)
    // Amplitude(s) for diagram number 9248
    FFV1_0( w_fp[662], w_fp[161], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9249 OF 15495 ***
    // Wavefunction(s) for diagram number 9249
    // (none)
    // Amplitude(s) for diagram number 9249
    FFV1_0( w_fp[662], w_fp[156], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];

    // *** DIAGRAM 9250 OF 15495 ***
    // Wavefunction(s) for diagram number 9250
    // (none)
    // Amplitude(s) for diagram number 9250
    FFV1_0( w_fp[245], w_fp[666], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9251 OF 15495 ***
    // Wavefunction(s) for diagram number 9251
    // (none)
    // Amplitude(s) for diagram number 9251
    FFV1_0( w_fp[39], w_fp[666], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9252 OF 15495 ***
    // Wavefunction(s) for diagram number 9252
    // (none)
    // Amplitude(s) for diagram number 9252
    FFV1_0( w_fp[157], w_fp[666], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];

    // *** DIAGRAM 9253 OF 15495 ***
    // Wavefunction(s) for diagram number 9253
    // (none)
    // Amplitude(s) for diagram number 9253
    VVV1_0( w_fp[522], w_fp[497], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9254 OF 15495 ***
    // Wavefunction(s) for diagram number 9254
    // (none)
    // Amplitude(s) for diagram number 9254
    FFV1_0( w_fp[39], w_fp[156], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];

    // *** DIAGRAM 9255 OF 15495 ***
    // Wavefunction(s) for diagram number 9255
    // (none)
    // Amplitude(s) for diagram number 9255
    FFV1_0( w_fp[157], w_fp[161], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];

    // *** DIAGRAM 9256 OF 15495 ***
    // Wavefunction(s) for diagram number 9256
    // (none)
    // Amplitude(s) for diagram number 9256
    VVV1_0( w_fp[0], w_fp[497], w_fp[125], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9257 OF 15495 ***
    // Wavefunction(s) for diagram number 9257
    // (none)
    // Amplitude(s) for diagram number 9257
    FFV1_0( w_fp[245], w_fp[161], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9258 OF 15495 ***
    // Wavefunction(s) for diagram number 9258
    // (none)
    // Amplitude(s) for diagram number 9258
    FFV1_0( w_fp[39], w_fp[187], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9259 OF 15495 ***
    // Wavefunction(s) for diagram number 9259
    // (none)
    // Amplitude(s) for diagram number 9259
    FFV1_0( w_fp[157], w_fp[156], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9260 OF 15495 ***
    // Wavefunction(s) for diagram number 9260
    // (none)
    // Amplitude(s) for diagram number 9260
    FFV1_0( w_fp[662], w_fp[158], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9261 OF 15495 ***
    // Wavefunction(s) for diagram number 9261
    // (none)
    // Amplitude(s) for diagram number 9261
    FFV1_0( w_fp[662], w_fp[189], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9262 OF 15495 ***
    // Wavefunction(s) for diagram number 9262
    // (none)
    // Amplitude(s) for diagram number 9262
    FFV1_0( w_fp[662], w_fp[156], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 9263 OF 15495 ***
    // Wavefunction(s) for diagram number 9263
    // (none)
    // Amplitude(s) for diagram number 9263
    FFV1_0( w_fp[29], w_fp[666], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9264 OF 15495 ***
    // Wavefunction(s) for diagram number 9264
    // (none)
    // Amplitude(s) for diagram number 9264
    FFV1_0( w_fp[246], w_fp[666], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9265 OF 15495 ***
    // Wavefunction(s) for diagram number 9265
    // (none)
    // Amplitude(s) for diagram number 9265
    FFV1_0( w_fp[157], w_fp[666], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];

    // *** DIAGRAM 9266 OF 15495 ***
    // Wavefunction(s) for diagram number 9266
    // (none)
    // Amplitude(s) for diagram number 9266
    VVV1_0( w_fp[516], w_fp[497], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9267 OF 15495 ***
    // Wavefunction(s) for diagram number 9267
    // (none)
    // Amplitude(s) for diagram number 9267
    FFV1_0( w_fp[29], w_fp[156], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];

    // *** DIAGRAM 9268 OF 15495 ***
    // Wavefunction(s) for diagram number 9268
    // (none)
    // Amplitude(s) for diagram number 9268
    FFV1_0( w_fp[157], w_fp[158], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];

    // *** DIAGRAM 9269 OF 15495 ***
    // Wavefunction(s) for diagram number 9269
    // (none)
    // Amplitude(s) for diagram number 9269
    VVV1_0( w_fp[0], w_fp[497], w_fp[131], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9270 OF 15495 ***
    // Wavefunction(s) for diagram number 9270
    // (none)
    // Amplitude(s) for diagram number 9270
    FFV1_0( w_fp[29], w_fp[189], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9271 OF 15495 ***
    // Wavefunction(s) for diagram number 9271
    // (none)
    // Amplitude(s) for diagram number 9271
    FFV1_0( w_fp[246], w_fp[158], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9272 OF 15495 ***
    // Wavefunction(s) for diagram number 9272
    // (none)
    // Amplitude(s) for diagram number 9272
    FFV1_0( w_fp[157], w_fp[156], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9273 OF 15495 ***
    // Wavefunction(s) for diagram number 9273
    // (none)
    // Amplitude(s) for diagram number 9273
    FFV1_0( w_fp[662], w_fp[156], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    FFV1_0( w_fp[662], w_fp[156], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    FFV1_0( w_fp[662], w_fp[156], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];

    // *** DIAGRAM 9274 OF 15495 ***
    // Wavefunction(s) for diagram number 9274
    // (none)
    // Amplitude(s) for diagram number 9274
    FFV1_0( w_fp[157], w_fp[666], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[666], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[666], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];

    // *** DIAGRAM 9275 OF 15495 ***
    // Wavefunction(s) for diagram number 9275
    // (none)
    // Amplitude(s) for diagram number 9275
    FFV1_0( w_fp[157], w_fp[156], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9276 OF 15495 ***
    // Wavefunction(s) for diagram number 9276
    // (none)
    // Amplitude(s) for diagram number 9276
    FFV1_0( w_fp[663], w_fp[190], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] -= amp_sv[0];

    // *** DIAGRAM 9277 OF 15495 ***
    // Wavefunction(s) for diagram number 9277
    // (none)
    // Amplitude(s) for diagram number 9277
    FFV1_0( w_fp[664], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] -= amp_sv[0];

    // *** DIAGRAM 9278 OF 15495 ***
    // Wavefunction(s) for diagram number 9278
    FFV1_2( w_fp[662], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[497] );
    // Amplitude(s) for diagram number 9278
    FFV1_0( w_fp[497], w_fp[191], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] -= amp_sv[0];

    // *** DIAGRAM 9279 OF 15495 ***
    // Wavefunction(s) for diagram number 9279
    // (none)
    // Amplitude(s) for diagram number 9279
    FFV1_0( w_fp[664], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] -= amp_sv[0];

    // *** DIAGRAM 9280 OF 15495 ***
    // Wavefunction(s) for diagram number 9280
    // (none)
    // Amplitude(s) for diagram number 9280
    FFV1_0( w_fp[497], w_fp[193], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 9281 OF 15495 ***
    // Wavefunction(s) for diagram number 9281
    // (none)
    // Amplitude(s) for diagram number 9281
    FFV1_0( w_fp[663], w_fp[193], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] -= amp_sv[0];

    // *** DIAGRAM 9282 OF 15495 ***
    // Wavefunction(s) for diagram number 9282
    FFV1_1( w_fp[169], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[676] );
    FFV1_1( w_fp[676], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[677] );
    // Amplitude(s) for diagram number 9282
    FFV1_0( w_fp[94], w_fp[677], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] -= amp_sv[0];

    // *** DIAGRAM 9283 OF 15495 ***
    // Wavefunction(s) for diagram number 9283
    FFV1_1( w_fp[676], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[678] );
    // Amplitude(s) for diagram number 9283
    FFV1_0( w_fp[94], w_fp[678], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[383] -= amp_sv[0];

    // *** DIAGRAM 9284 OF 15495 ***
    // Wavefunction(s) for diagram number 9284
    FFV1_1( w_fp[676], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[679] );
    // Amplitude(s) for diagram number 9284
    FFV1_0( w_fp[39], w_fp[679], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];

    // *** DIAGRAM 9285 OF 15495 ***
    // Wavefunction(s) for diagram number 9285
    // (none)
    // Amplitude(s) for diagram number 9285
    FFV1_0( w_fp[39], w_fp[678], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] -= amp_sv[0];

    // *** DIAGRAM 9286 OF 15495 ***
    // Wavefunction(s) for diagram number 9286
    // (none)
    // Amplitude(s) for diagram number 9286
    FFV1_0( w_fp[77], w_fp[679], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];

    // *** DIAGRAM 9287 OF 15495 ***
    // Wavefunction(s) for diagram number 9287
    // (none)
    // Amplitude(s) for diagram number 9287
    FFV1_0( w_fp[77], w_fp[677], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= amp_sv[0];

    // *** DIAGRAM 9288 OF 15495 ***
    // Wavefunction(s) for diagram number 9288
    FFV1_2( w_fp[94], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[680] );
    // Amplitude(s) for diagram number 9288
    FFV1_0( w_fp[680], w_fp[191], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] -= amp_sv[0];

    // *** DIAGRAM 9289 OF 15495 ***
    // Wavefunction(s) for diagram number 9289
    FFV1_1( w_fp[191], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[681] );
    // Amplitude(s) for diagram number 9289
    FFV1_0( w_fp[94], w_fp[681], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[437] -= amp_sv[0];

    // *** DIAGRAM 9290 OF 15495 ***
    // Wavefunction(s) for diagram number 9290
    // (none)
    // Amplitude(s) for diagram number 9290
    FFV1_0( w_fp[680], w_fp[193], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 9291 OF 15495 ***
    // Wavefunction(s) for diagram number 9291
    FFV1_1( w_fp[193], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[682] );
    // Amplitude(s) for diagram number 9291
    FFV1_0( w_fp[94], w_fp[682], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] -= amp_sv[0];

    // *** DIAGRAM 9292 OF 15495 ***
    // Wavefunction(s) for diagram number 9292
    // (none)
    // Amplitude(s) for diagram number 9292
    FFV1_0( w_fp[673], w_fp[190], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] -= amp_sv[0];

    // *** DIAGRAM 9293 OF 15495 ***
    // Wavefunction(s) for diagram number 9293
    FFV1_1( w_fp[190], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[683] );
    // Amplitude(s) for diagram number 9293
    FFV1_0( w_fp[39], w_fp[683], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[413] -= amp_sv[0];

    // *** DIAGRAM 9294 OF 15495 ***
    // Wavefunction(s) for diagram number 9294
    // (none)
    // Amplitude(s) for diagram number 9294
    FFV1_0( w_fp[673], w_fp[193], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[469] -= amp_sv[0];

    // *** DIAGRAM 9295 OF 15495 ***
    // Wavefunction(s) for diagram number 9295
    // (none)
    // Amplitude(s) for diagram number 9295
    FFV1_0( w_fp[39], w_fp[682], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] -= amp_sv[0];

    // *** DIAGRAM 9296 OF 15495 ***
    // Wavefunction(s) for diagram number 9296
    // (none)
    // Amplitude(s) for diagram number 9296
    FFV1_0( w_fp[675], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[421] -= amp_sv[0];

    // *** DIAGRAM 9297 OF 15495 ***
    // Wavefunction(s) for diagram number 9297
    // (none)
    // Amplitude(s) for diagram number 9297
    FFV1_0( w_fp[77], w_fp[683], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] -= amp_sv[0];

    // *** DIAGRAM 9298 OF 15495 ***
    // Wavefunction(s) for diagram number 9298
    // (none)
    // Amplitude(s) for diagram number 9298
    FFV1_0( w_fp[675], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[445] -= amp_sv[0];

    // *** DIAGRAM 9299 OF 15495 ***
    // Wavefunction(s) for diagram number 9299
    // (none)
    // Amplitude(s) for diagram number 9299
    FFV1_0( w_fp[77], w_fp[681], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] -= amp_sv[0];

    // *** DIAGRAM 9300 OF 15495 ***
    // Wavefunction(s) for diagram number 9300
    // (none)
    // Amplitude(s) for diagram number 9300
    FFV1_0( w_fp[662], w_fp[207], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9301 OF 15495 ***
    // Wavefunction(s) for diagram number 9301
    // (none)
    // Amplitude(s) for diagram number 9301
    FFV1_0( w_fp[662], w_fp[193], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9302 OF 15495 ***
    // Wavefunction(s) for diagram number 9302
    // (none)
    // Amplitude(s) for diagram number 9302
    FFV1_0( w_fp[662], w_fp[169], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 9303 OF 15495 ***
    // Wavefunction(s) for diagram number 9303
    // (none)
    // Amplitude(s) for diagram number 9303
    FFV1_0( w_fp[110], w_fp[676], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9304 OF 15495 ***
    // Wavefunction(s) for diagram number 9304
    // (none)
    // Amplitude(s) for diagram number 9304
    FFV1_0( w_fp[77], w_fp[676], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9305 OF 15495 ***
    // Wavefunction(s) for diagram number 9305
    // (none)
    // Amplitude(s) for diagram number 9305
    FFV1_0( w_fp[157], w_fp[676], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];

    // *** DIAGRAM 9306 OF 15495 ***
    // Wavefunction(s) for diagram number 9306
    // (none)
    // Amplitude(s) for diagram number 9306
    VVV1_0( w_fp[576], w_fp[540], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9307 OF 15495 ***
    // Wavefunction(s) for diagram number 9307
    // (none)
    // Amplitude(s) for diagram number 9307
    FFV1_0( w_fp[77], w_fp[169], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];

    // *** DIAGRAM 9308 OF 15495 ***
    // Wavefunction(s) for diagram number 9308
    // (none)
    // Amplitude(s) for diagram number 9308
    FFV1_0( w_fp[157], w_fp[193], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 9309 OF 15495 ***
    // Wavefunction(s) for diagram number 9309
    // (none)
    // Amplitude(s) for diagram number 9309
    VVV1_0( w_fp[0], w_fp[540], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9310 OF 15495 ***
    // Wavefunction(s) for diagram number 9310
    // (none)
    // Amplitude(s) for diagram number 9310
    FFV1_0( w_fp[110], w_fp[193], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9311 OF 15495 ***
    // Wavefunction(s) for diagram number 9311
    // (none)
    // Amplitude(s) for diagram number 9311
    FFV1_0( w_fp[77], w_fp[207], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9312 OF 15495 ***
    // Wavefunction(s) for diagram number 9312
    // (none)
    // Amplitude(s) for diagram number 9312
    FFV1_0( w_fp[157], w_fp[169], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9313 OF 15495 ***
    // Wavefunction(s) for diagram number 9313
    // (none)
    // Amplitude(s) for diagram number 9313
    FFV1_0( w_fp[662], w_fp[209], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9314 OF 15495 ***
    // Wavefunction(s) for diagram number 9314
    // (none)
    // Amplitude(s) for diagram number 9314
    FFV1_0( w_fp[662], w_fp[191], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9315 OF 15495 ***
    // Wavefunction(s) for diagram number 9315
    // (none)
    // Amplitude(s) for diagram number 9315
    FFV1_0( w_fp[662], w_fp[169], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];

    // *** DIAGRAM 9316 OF 15495 ***
    // Wavefunction(s) for diagram number 9316
    // (none)
    // Amplitude(s) for diagram number 9316
    FFV1_0( w_fp[120], w_fp[676], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9317 OF 15495 ***
    // Wavefunction(s) for diagram number 9317
    // (none)
    // Amplitude(s) for diagram number 9317
    FFV1_0( w_fp[39], w_fp[676], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9318 OF 15495 ***
    // Wavefunction(s) for diagram number 9318
    // (none)
    // Amplitude(s) for diagram number 9318
    FFV1_0( w_fp[157], w_fp[676], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];

    // *** DIAGRAM 9319 OF 15495 ***
    // Wavefunction(s) for diagram number 9319
    // (none)
    // Amplitude(s) for diagram number 9319
    VVV1_0( w_fp[518], w_fp[540], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9320 OF 15495 ***
    // Wavefunction(s) for diagram number 9320
    // (none)
    // Amplitude(s) for diagram number 9320
    FFV1_0( w_fp[39], w_fp[169], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];

    // *** DIAGRAM 9321 OF 15495 ***
    // Wavefunction(s) for diagram number 9321
    // (none)
    // Amplitude(s) for diagram number 9321
    FFV1_0( w_fp[157], w_fp[191], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];

    // *** DIAGRAM 9322 OF 15495 ***
    // Wavefunction(s) for diagram number 9322
    // (none)
    // Amplitude(s) for diagram number 9322
    VVV1_0( w_fp[0], w_fp[540], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9323 OF 15495 ***
    // Wavefunction(s) for diagram number 9323
    // (none)
    // Amplitude(s) for diagram number 9323
    FFV1_0( w_fp[120], w_fp[191], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9324 OF 15495 ***
    // Wavefunction(s) for diagram number 9324
    // (none)
    // Amplitude(s) for diagram number 9324
    FFV1_0( w_fp[39], w_fp[209], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9325 OF 15495 ***
    // Wavefunction(s) for diagram number 9325
    // (none)
    // Amplitude(s) for diagram number 9325
    FFV1_0( w_fp[157], w_fp[169], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[529], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9326 OF 15495 ***
    // Wavefunction(s) for diagram number 9326
    // (none)
    // Amplitude(s) for diagram number 9326
    FFV1_0( w_fp[662], w_fp[190], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9327 OF 15495 ***
    // Wavefunction(s) for diagram number 9327
    // (none)
    // Amplitude(s) for diagram number 9327
    FFV1_0( w_fp[662], w_fp[210], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9328 OF 15495 ***
    // Wavefunction(s) for diagram number 9328
    // (none)
    // Amplitude(s) for diagram number 9328
    FFV1_0( w_fp[662], w_fp[169], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 9329 OF 15495 ***
    // Wavefunction(s) for diagram number 9329
    // (none)
    // Amplitude(s) for diagram number 9329
    FFV1_0( w_fp[94], w_fp[676], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9330 OF 15495 ***
    // Wavefunction(s) for diagram number 9330
    // (none)
    // Amplitude(s) for diagram number 9330
    FFV1_0( w_fp[246], w_fp[676], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9331 OF 15495 ***
    // Wavefunction(s) for diagram number 9331
    // (none)
    // Amplitude(s) for diagram number 9331
    FFV1_0( w_fp[157], w_fp[676], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];

    // *** DIAGRAM 9332 OF 15495 ***
    // Wavefunction(s) for diagram number 9332
    // (none)
    // Amplitude(s) for diagram number 9332
    VVV1_0( w_fp[516], w_fp[540], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9333 OF 15495 ***
    // Wavefunction(s) for diagram number 9333
    // (none)
    // Amplitude(s) for diagram number 9333
    FFV1_0( w_fp[94], w_fp[169], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 9334 OF 15495 ***
    // Wavefunction(s) for diagram number 9334
    // (none)
    // Amplitude(s) for diagram number 9334
    FFV1_0( w_fp[157], w_fp[190], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];

    // *** DIAGRAM 9335 OF 15495 ***
    // Wavefunction(s) for diagram number 9335
    // (none)
    // Amplitude(s) for diagram number 9335
    VVV1_0( w_fp[0], w_fp[540], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9336 OF 15495 ***
    // Wavefunction(s) for diagram number 9336
    // (none)
    // Amplitude(s) for diagram number 9336
    FFV1_0( w_fp[94], w_fp[210], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9337 OF 15495 ***
    // Wavefunction(s) for diagram number 9337
    // (none)
    // Amplitude(s) for diagram number 9337
    FFV1_0( w_fp[246], w_fp[190], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9338 OF 15495 ***
    // Wavefunction(s) for diagram number 9338
    // (none)
    // Amplitude(s) for diagram number 9338
    FFV1_0( w_fp[157], w_fp[169], w_fp[468], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9339 OF 15495 ***
    // Wavefunction(s) for diagram number 9339
    // (none)
    // Amplitude(s) for diagram number 9339
    FFV1_0( w_fp[662], w_fp[169], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    FFV1_0( w_fp[662], w_fp[169], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    FFV1_0( w_fp[662], w_fp[169], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 9340 OF 15495 ***
    // Wavefunction(s) for diagram number 9340
    // (none)
    // Amplitude(s) for diagram number 9340
    FFV1_0( w_fp[157], w_fp[676], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[676], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[676], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];

    // *** DIAGRAM 9341 OF 15495 ***
    // Wavefunction(s) for diagram number 9341
    // (none)
    // Amplitude(s) for diagram number 9341
    FFV1_0( w_fp[157], w_fp[169], w_fp[593], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9342 OF 15495 ***
    // Wavefunction(s) for diagram number 9342
    // (none)
    // Amplitude(s) for diagram number 9342
    FFV1_0( w_fp[665], w_fp[211], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] -= amp_sv[0];

    // *** DIAGRAM 9343 OF 15495 ***
    // Wavefunction(s) for diagram number 9343
    // (none)
    // Amplitude(s) for diagram number 9343
    FFV1_0( w_fp[664], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] -= amp_sv[0];

    // *** DIAGRAM 9344 OF 15495 ***
    // Wavefunction(s) for diagram number 9344
    // (none)
    // Amplitude(s) for diagram number 9344
    FFV1_0( w_fp[497], w_fp[212], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[574] -= amp_sv[0];

    // *** DIAGRAM 9345 OF 15495 ***
    // Wavefunction(s) for diagram number 9345
    // (none)
    // Amplitude(s) for diagram number 9345
    FFV1_0( w_fp[664], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[568] -= amp_sv[0];

    // *** DIAGRAM 9346 OF 15495 ***
    // Wavefunction(s) for diagram number 9346
    // (none)
    // Amplitude(s) for diagram number 9346
    FFV1_0( w_fp[497], w_fp[213], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 9347 OF 15495 ***
    // Wavefunction(s) for diagram number 9347
    // (none)
    // Amplitude(s) for diagram number 9347
    FFV1_0( w_fp[665], w_fp[213], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] -= amp_sv[0];

    // *** DIAGRAM 9348 OF 15495 ***
    // Wavefunction(s) for diagram number 9348
    FFV1_1( w_fp[197], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[540] );
    FFV1_1( w_fp[540], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[684] );
    // Amplitude(s) for diagram number 9348
    FFV1_0( w_fp[94], w_fp[684], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] -= amp_sv[0];

    // *** DIAGRAM 9349 OF 15495 ***
    // Wavefunction(s) for diagram number 9349
    FFV1_1( w_fp[540], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[685] );
    // Amplitude(s) for diagram number 9349
    FFV1_0( w_fp[94], w_fp[685], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[503] -= amp_sv[0];

    // *** DIAGRAM 9350 OF 15495 ***
    // Wavefunction(s) for diagram number 9350
    FFV1_1( w_fp[540], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[686] );
    // Amplitude(s) for diagram number 9350
    FFV1_0( w_fp[29], w_fp[686], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= amp_sv[0];

    // *** DIAGRAM 9351 OF 15495 ***
    // Wavefunction(s) for diagram number 9351
    // (none)
    // Amplitude(s) for diagram number 9351
    FFV1_0( w_fp[29], w_fp[685], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] -= amp_sv[0];

    // *** DIAGRAM 9352 OF 15495 ***
    // Wavefunction(s) for diagram number 9352
    // (none)
    // Amplitude(s) for diagram number 9352
    FFV1_0( w_fp[77], w_fp[686], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];

    // *** DIAGRAM 9353 OF 15495 ***
    // Wavefunction(s) for diagram number 9353
    // (none)
    // Amplitude(s) for diagram number 9353
    FFV1_0( w_fp[77], w_fp[684], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] -= amp_sv[0];

    // *** DIAGRAM 9354 OF 15495 ***
    // Wavefunction(s) for diagram number 9354
    // (none)
    // Amplitude(s) for diagram number 9354
    FFV1_0( w_fp[680], w_fp[212], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[571] -= amp_sv[0];

    // *** DIAGRAM 9355 OF 15495 ***
    // Wavefunction(s) for diagram number 9355
    FFV1_1( w_fp[212], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[687] );
    // Amplitude(s) for diagram number 9355
    FFV1_0( w_fp[94], w_fp[687], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[557] -= amp_sv[0];

    // *** DIAGRAM 9356 OF 15495 ***
    // Wavefunction(s) for diagram number 9356
    // (none)
    // Amplitude(s) for diagram number 9356
    FFV1_0( w_fp[680], w_fp[213], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[595] -= amp_sv[0];

    // *** DIAGRAM 9357 OF 15495 ***
    // Wavefunction(s) for diagram number 9357
    FFV1_1( w_fp[213], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[688] );
    // Amplitude(s) for diagram number 9357
    FFV1_0( w_fp[94], w_fp[688], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] -= amp_sv[0];

    // *** DIAGRAM 9358 OF 15495 ***
    // Wavefunction(s) for diagram number 9358
    // (none)
    // Amplitude(s) for diagram number 9358
    FFV1_0( w_fp[670], w_fp[211], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[547] -= amp_sv[0];

    // *** DIAGRAM 9359 OF 15495 ***
    // Wavefunction(s) for diagram number 9359
    FFV1_1( w_fp[211], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[689] );
    // Amplitude(s) for diagram number 9359
    FFV1_0( w_fp[29], w_fp[689], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[533] -= amp_sv[0];

    // *** DIAGRAM 9360 OF 15495 ***
    // Wavefunction(s) for diagram number 9360
    // (none)
    // Amplitude(s) for diagram number 9360
    FFV1_0( w_fp[670], w_fp[213], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[589] -= amp_sv[0];

    // *** DIAGRAM 9361 OF 15495 ***
    // Wavefunction(s) for diagram number 9361
    // (none)
    // Amplitude(s) for diagram number 9361
    FFV1_0( w_fp[29], w_fp[688], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= amp_sv[0];

    // *** DIAGRAM 9362 OF 15495 ***
    // Wavefunction(s) for diagram number 9362
    // (none)
    // Amplitude(s) for diagram number 9362
    FFV1_0( w_fp[675], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[541] -= amp_sv[0];

    // *** DIAGRAM 9363 OF 15495 ***
    // Wavefunction(s) for diagram number 9363
    // (none)
    // Amplitude(s) for diagram number 9363
    FFV1_0( w_fp[77], w_fp[689], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[531] -= amp_sv[0];

    // *** DIAGRAM 9364 OF 15495 ***
    // Wavefunction(s) for diagram number 9364
    // (none)
    // Amplitude(s) for diagram number 9364
    FFV1_0( w_fp[675], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[565] -= amp_sv[0];

    // *** DIAGRAM 9365 OF 15495 ***
    // Wavefunction(s) for diagram number 9365
    // (none)
    // Amplitude(s) for diagram number 9365
    FFV1_0( w_fp[77], w_fp[687], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[555] -= amp_sv[0];

    // *** DIAGRAM 9366 OF 15495 ***
    // Wavefunction(s) for diagram number 9366
    // (none)
    // Amplitude(s) for diagram number 9366
    FFV1_0( w_fp[662], w_fp[222], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9367 OF 15495 ***
    // Wavefunction(s) for diagram number 9367
    // (none)
    // Amplitude(s) for diagram number 9367
    FFV1_0( w_fp[662], w_fp[213], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9368 OF 15495 ***
    // Wavefunction(s) for diagram number 9368
    // (none)
    // Amplitude(s) for diagram number 9368
    FFV1_0( w_fp[662], w_fp[197], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 9369 OF 15495 ***
    // Wavefunction(s) for diagram number 9369
    // (none)
    // Amplitude(s) for diagram number 9369
    FFV1_0( w_fp[92], w_fp[540], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9370 OF 15495 ***
    // Wavefunction(s) for diagram number 9370
    // (none)
    // Amplitude(s) for diagram number 9370
    FFV1_0( w_fp[77], w_fp[540], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9371 OF 15495 ***
    // Wavefunction(s) for diagram number 9371
    // (none)
    // Amplitude(s) for diagram number 9371
    FFV1_0( w_fp[157], w_fp[540], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];

    // *** DIAGRAM 9372 OF 15495 ***
    // Wavefunction(s) for diagram number 9372
    // (none)
    // Amplitude(s) for diagram number 9372
    VVV1_0( w_fp[514], w_fp[542], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9373 OF 15495 ***
    // Wavefunction(s) for diagram number 9373
    // (none)
    // Amplitude(s) for diagram number 9373
    FFV1_0( w_fp[77], w_fp[197], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];

    // *** DIAGRAM 9374 OF 15495 ***
    // Wavefunction(s) for diagram number 9374
    // (none)
    // Amplitude(s) for diagram number 9374
    FFV1_0( w_fp[157], w_fp[213], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 9375 OF 15495 ***
    // Wavefunction(s) for diagram number 9375
    // (none)
    // Amplitude(s) for diagram number 9375
    VVV1_0( w_fp[0], w_fp[542], w_fp[69], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9376 OF 15495 ***
    // Wavefunction(s) for diagram number 9376
    // (none)
    // Amplitude(s) for diagram number 9376
    FFV1_0( w_fp[92], w_fp[213], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9377 OF 15495 ***
    // Wavefunction(s) for diagram number 9377
    // (none)
    // Amplitude(s) for diagram number 9377
    FFV1_0( w_fp[77], w_fp[222], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9378 OF 15495 ***
    // Wavefunction(s) for diagram number 9378
    // (none)
    // Amplitude(s) for diagram number 9378
    FFV1_0( w_fp[157], w_fp[197], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9379 OF 15495 ***
    // Wavefunction(s) for diagram number 9379
    // (none)
    // Amplitude(s) for diagram number 9379
    FFV1_0( w_fp[662], w_fp[223], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9380 OF 15495 ***
    // Wavefunction(s) for diagram number 9380
    // (none)
    // Amplitude(s) for diagram number 9380
    FFV1_0( w_fp[662], w_fp[212], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9381 OF 15495 ***
    // Wavefunction(s) for diagram number 9381
    // (none)
    // Amplitude(s) for diagram number 9381
    FFV1_0( w_fp[662], w_fp[197], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];

    // *** DIAGRAM 9382 OF 15495 ***
    // Wavefunction(s) for diagram number 9382
    // (none)
    // Amplitude(s) for diagram number 9382
    FFV1_0( w_fp[120], w_fp[540], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9383 OF 15495 ***
    // Wavefunction(s) for diagram number 9383
    // (none)
    // Amplitude(s) for diagram number 9383
    FFV1_0( w_fp[29], w_fp[540], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9384 OF 15495 ***
    // Wavefunction(s) for diagram number 9384
    // (none)
    // Amplitude(s) for diagram number 9384
    FFV1_0( w_fp[157], w_fp[540], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];

    // *** DIAGRAM 9385 OF 15495 ***
    // Wavefunction(s) for diagram number 9385
    // (none)
    // Amplitude(s) for diagram number 9385
    VVV1_0( w_fp[518], w_fp[542], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9386 OF 15495 ***
    // Wavefunction(s) for diagram number 9386
    // (none)
    // Amplitude(s) for diagram number 9386
    FFV1_0( w_fp[29], w_fp[197], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];

    // *** DIAGRAM 9387 OF 15495 ***
    // Wavefunction(s) for diagram number 9387
    // (none)
    // Amplitude(s) for diagram number 9387
    FFV1_0( w_fp[157], w_fp[212], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];

    // *** DIAGRAM 9388 OF 15495 ***
    // Wavefunction(s) for diagram number 9388
    // (none)
    // Amplitude(s) for diagram number 9388
    VVV1_0( w_fp[0], w_fp[542], w_fp[103], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9389 OF 15495 ***
    // Wavefunction(s) for diagram number 9389
    // (none)
    // Amplitude(s) for diagram number 9389
    FFV1_0( w_fp[120], w_fp[212], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9390 OF 15495 ***
    // Wavefunction(s) for diagram number 9390
    // (none)
    // Amplitude(s) for diagram number 9390
    FFV1_0( w_fp[29], w_fp[223], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9391 OF 15495 ***
    // Wavefunction(s) for diagram number 9391
    // (none)
    // Amplitude(s) for diagram number 9391
    FFV1_0( w_fp[157], w_fp[197], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[471], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9392 OF 15495 ***
    // Wavefunction(s) for diagram number 9392
    // (none)
    // Amplitude(s) for diagram number 9392
    FFV1_0( w_fp[662], w_fp[211], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9393 OF 15495 ***
    // Wavefunction(s) for diagram number 9393
    // (none)
    // Amplitude(s) for diagram number 9393
    FFV1_0( w_fp[662], w_fp[224], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9394 OF 15495 ***
    // Wavefunction(s) for diagram number 9394
    // (none)
    // Amplitude(s) for diagram number 9394
    FFV1_0( w_fp[662], w_fp[197], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 9395 OF 15495 ***
    // Wavefunction(s) for diagram number 9395
    // (none)
    // Amplitude(s) for diagram number 9395
    FFV1_0( w_fp[94], w_fp[540], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9396 OF 15495 ***
    // Wavefunction(s) for diagram number 9396
    // (none)
    // Amplitude(s) for diagram number 9396
    FFV1_0( w_fp[245], w_fp[540], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9397 OF 15495 ***
    // Wavefunction(s) for diagram number 9397
    // (none)
    // Amplitude(s) for diagram number 9397
    FFV1_0( w_fp[157], w_fp[540], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];

    // *** DIAGRAM 9398 OF 15495 ***
    // Wavefunction(s) for diagram number 9398
    // (none)
    // Amplitude(s) for diagram number 9398
    VVV1_0( w_fp[522], w_fp[542], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 9399 OF 15495 ***
    // Wavefunction(s) for diagram number 9399
    // (none)
    // Amplitude(s) for diagram number 9399
    FFV1_0( w_fp[94], w_fp[197], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];

    // *** DIAGRAM 9400 OF 15495 ***
    // Wavefunction(s) for diagram number 9400
    // (none)
    // Amplitude(s) for diagram number 9400
    FFV1_0( w_fp[157], w_fp[211], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 449 );
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 540 );
    storeWf( wfs, w_cx, nevt, 645 );
    storeWf( wfs, w_cx, nevt, 646 );
    storeWf( wfs, w_cx, nevt, 647 );
    storeWf( wfs, w_cx, nevt, 648 );
    storeWf( wfs, w_cx, nevt, 649 );
    storeWf( wfs, w_cx, nevt, 650 );
    storeWf( wfs, w_cx, nevt, 651 );
    storeWf( wfs, w_cx, nevt, 652 );
    storeWf( wfs, w_cx, nevt, 653 );
    storeWf( wfs, w_cx, nevt, 654 );
    storeWf( wfs, w_cx, nevt, 655 );
    storeWf( wfs, w_cx, nevt, 656 );
    storeWf( wfs, w_cx, nevt, 657 );
    storeWf( wfs, w_cx, nevt, 658 );
    storeWf( wfs, w_cx, nevt, 659 );
    storeWf( wfs, w_cx, nevt, 660 );
    storeWf( wfs, w_cx, nevt, 661 );
    storeWf( wfs, w_cx, nevt, 662 );
    storeWf( wfs, w_cx, nevt, 663 );
    storeWf( wfs, w_cx, nevt, 664 );
    storeWf( wfs, w_cx, nevt, 665 );
    storeWf( wfs, w_cx, nevt, 666 );
    storeWf( wfs, w_cx, nevt, 667 );
    storeWf( wfs, w_cx, nevt, 668 );
    storeWf( wfs, w_cx, nevt, 669 );
    storeWf( wfs, w_cx, nevt, 670 );
    storeWf( wfs, w_cx, nevt, 671 );
    storeWf( wfs, w_cx, nevt, 672 );
    storeWf( wfs, w_cx, nevt, 673 );
    storeWf( wfs, w_cx, nevt, 674 );
    storeWf( wfs, w_cx, nevt, 675 );
    storeWf( wfs, w_cx, nevt, 676 );
    storeWf( wfs, w_cx, nevt, 677 );
    storeWf( wfs, w_cx, nevt, 678 );
    storeWf( wfs, w_cx, nevt, 679 );
    storeWf( wfs, w_cx, nevt, 680 );
    storeWf( wfs, w_cx, nevt, 681 );
    storeWf( wfs, w_cx, nevt, 682 );
    storeWf( wfs, w_cx, nevt, 683 );
    storeWf( wfs, w_cx, nevt, 684 );
    storeWf( wfs, w_cx, nevt, 685 );
    storeWf( wfs, w_cx, nevt, 686 );
    storeWf( wfs, w_cx, nevt, 687 );
    storeWf( wfs, w_cx, nevt, 688 );
    storeWf( wfs, w_cx, nevt, 689 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
