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
  diagramgroup59( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 9 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
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
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 273 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 281 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 328 );
    retrieveWf( wfs, w_cx, nevt, 329 );
    retrieveWf( wfs, w_cx, nevt, 332 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 335 );
    retrieveWf( wfs, w_cx, nevt, 336 );
    retrieveWf( wfs, w_cx, nevt, 337 );
    retrieveWf( wfs, w_cx, nevt, 338 );
    retrieveWf( wfs, w_cx, nevt, 339 );
    retrieveWf( wfs, w_cx, nevt, 340 );
    retrieveWf( wfs, w_cx, nevt, 341 );
    retrieveWf( wfs, w_cx, nevt, 342 );
    retrieveWf( wfs, w_cx, nevt, 343 );
    retrieveWf( wfs, w_cx, nevt, 344 );
    retrieveWf( wfs, w_cx, nevt, 346 );
    retrieveWf( wfs, w_cx, nevt, 347 );
    retrieveWf( wfs, w_cx, nevt, 348 );
    retrieveWf( wfs, w_cx, nevt, 350 );
    retrieveWf( wfs, w_cx, nevt, 351 );
    retrieveWf( wfs, w_cx, nevt, 489 );
    retrieveWf( wfs, w_cx, nevt, 498 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 503 );
    retrieveWf( wfs, w_cx, nevt, 504 );
    retrieveWf( wfs, w_cx, nevt, 506 );
    retrieveWf( wfs, w_cx, nevt, 507 );
    retrieveWf( wfs, w_cx, nevt, 508 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 548 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 589 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 665 );
    retrieveWf( wfs, w_cx, nevt, 666 );
    retrieveWf( wfs, w_cx, nevt, 667 );
    retrieveWf( wfs, w_cx, nevt, 669 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 671 );
    retrieveWf( wfs, w_cx, nevt, 674 );
    retrieveWf( wfs, w_cx, nevt, 676 );
    retrieveWf( wfs, w_cx, nevt, 677 );
    retrieveWf( wfs, w_cx, nevt, 679 );
    retrieveWf( wfs, w_cx, nevt, 681 );
    retrieveWf( wfs, w_cx, nevt, 683 );
    retrieveWf( wfs, w_cx, nevt, 684 );
    retrieveWf( wfs, w_cx, nevt, 686 );
    retrieveWf( wfs, w_cx, nevt, 687 );
    retrieveWf( wfs, w_cx, nevt, 689 );
    retrieveWf( wfs, w_cx, nevt, 703 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 707 );
    retrieveWf( wfs, w_cx, nevt, 708 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 742 );
    retrieveWf( wfs, w_cx, nevt, 743 );
    retrieveWf( wfs, w_cx, nevt, 744 );
    retrieveWf( wfs, w_cx, nevt, 745 );
    retrieveWf( wfs, w_cx, nevt, 748 );
#endif
#endif

    // *** DIAGRAM 11601 OF 15495 ***
    // Wavefunction(s) for diagram number 11601
    // (none)
    // Amplitude(s) for diagram number 11601
    VVV1_0( w_fp[325], w_fp[133], w_fp[665], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[134], w_fp[665], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[135], w_fp[665], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11602 OF 15495 ***
    // Wavefunction(s) for diagram number 11602
    // (none)
    // Amplitude(s) for diagram number 11602
    VVV1_0( w_fp[325], w_fp[9], w_fp[515], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[9], w_fp[589], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[9], w_fp[510], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11603 OF 15495 ***
    // Wavefunction(s) for diagram number 11603
    // (none)
    // Amplitude(s) for diagram number 11603
    VVVV1_0( w_fp[662], w_fp[164], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[164], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[164], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];

    // *** DIAGRAM 11604 OF 15495 ***
    // Wavefunction(s) for diagram number 11604
    // (none)
    // Amplitude(s) for diagram number 11604
    VVV1_0( w_fp[164], w_fp[6], w_fp[706], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];

    // *** DIAGRAM 11605 OF 15495 ***
    // Wavefunction(s) for diagram number 11605
    // (none)
    // Amplitude(s) for diagram number 11605
    VVV1_0( w_fp[164], w_fp[5], w_fp[500], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];

    // *** DIAGRAM 11606 OF 15495 ***
    // Wavefunction(s) for diagram number 11606
    FFV1_2( w_fp[3], w_fp[662], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[554] );
    // Amplitude(s) for diagram number 11606
    FFV1_0( w_fp[554], w_fp[158], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];

    // *** DIAGRAM 11607 OF 15495 ***
    // Wavefunction(s) for diagram number 11607
    // (none)
    // Amplitude(s) for diagram number 11607
    FFV1_0( w_fp[3], w_fp[158], w_fp[500], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11608 OF 15495 ***
    // Wavefunction(s) for diagram number 11608
    // (none)
    // Amplitude(s) for diagram number 11608
    FFV1_0( w_fp[554], w_fp[161], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];

    // *** DIAGRAM 11609 OF 15495 ***
    // Wavefunction(s) for diagram number 11609
    // (none)
    // Amplitude(s) for diagram number 11609
    FFV1_0( w_fp[3], w_fp[161], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11610 OF 15495 ***
    // Wavefunction(s) for diagram number 11610
    // (none)
    // Amplitude(s) for diagram number 11610
    FFV1_0( w_fp[344], w_fp[669], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11611 OF 15495 ***
    // Wavefunction(s) for diagram number 11611
    // (none)
    // Amplitude(s) for diagram number 11611
    FFV1_0( w_fp[344], w_fp[667], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11612 OF 15495 ***
    // Wavefunction(s) for diagram number 11612
    // (none)
    // Amplitude(s) for diagram number 11612
    VVV1_0( w_fp[333], w_fp[6], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11613 OF 15495 ***
    // Wavefunction(s) for diagram number 11613
    // (none)
    // Amplitude(s) for diagram number 11613
    FFV1_0( w_fp[3], w_fp[667], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[253] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];

    // *** DIAGRAM 11614 OF 15495 ***
    // Wavefunction(s) for diagram number 11614
    // (none)
    // Amplitude(s) for diagram number 11614
    VVV1_0( w_fp[334], w_fp[5], w_fp[281], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11615 OF 15495 ***
    // Wavefunction(s) for diagram number 11615
    // (none)
    // Amplitude(s) for diagram number 11615
    FFV1_0( w_fp[3], w_fp[669], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[247] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];

    // *** DIAGRAM 11616 OF 15495 ***
    // Wavefunction(s) for diagram number 11616
    // (none)
    // Amplitude(s) for diagram number 11616
    FFV1_0( w_fp[3], w_fp[666], w_fp[341], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[342], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[666], w_fp[343], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11617 OF 15495 ***
    // Wavefunction(s) for diagram number 11617
    FFV1_2( w_fp[344], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[317] );
    // Amplitude(s) for diagram number 11617
    FFV1_0( w_fp[317], w_fp[158], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11618 OF 15495 ***
    // Wavefunction(s) for diagram number 11618
    // (none)
    // Amplitude(s) for diagram number 11618
    FFV1_0( w_fp[344], w_fp[674], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11619 OF 15495 ***
    // Wavefunction(s) for diagram number 11619
    // (none)
    // Amplitude(s) for diagram number 11619
    FFV1_0( w_fp[317], w_fp[161], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11620 OF 15495 ***
    // Wavefunction(s) for diagram number 11620
    // (none)
    // Amplitude(s) for diagram number 11620
    FFV1_0( w_fp[344], w_fp[671], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11621 OF 15495 ***
    // Wavefunction(s) for diagram number 11621
    // (none)
    // Amplitude(s) for diagram number 11621
    VVVV1_0( w_fp[0], w_fp[333], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[333], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[333], w_fp[164], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];

    // *** DIAGRAM 11622 OF 15495 ***
    // Wavefunction(s) for diagram number 11622
    // (none)
    // Amplitude(s) for diagram number 11622
    VVV1_0( w_fp[164], w_fp[6], w_fp[270], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];

    // *** DIAGRAM 11623 OF 15495 ***
    // Wavefunction(s) for diagram number 11623
    // (none)
    // Amplitude(s) for diagram number 11623
    VVV1_0( w_fp[333], w_fp[6], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];

    // *** DIAGRAM 11624 OF 15495 ***
    // Wavefunction(s) for diagram number 11624
    // (none)
    // Amplitude(s) for diagram number 11624
    FFV1_0( w_fp[3], w_fp[161], w_fp[270], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11625 OF 15495 ***
    // Wavefunction(s) for diagram number 11625
    // (none)
    // Amplitude(s) for diagram number 11625
    FFV1_0( w_fp[3], w_fp[671], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];

    // *** DIAGRAM 11626 OF 15495 ***
    // Wavefunction(s) for diagram number 11626
    // (none)
    // Amplitude(s) for diagram number 11626
    VVVV1_0( w_fp[0], w_fp[334], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[334], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[334], w_fp[164], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];

    // *** DIAGRAM 11627 OF 15495 ***
    // Wavefunction(s) for diagram number 11627
    // (none)
    // Amplitude(s) for diagram number 11627
    VVV1_0( w_fp[164], w_fp[5], w_fp[503], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];

    // *** DIAGRAM 11628 OF 15495 ***
    // Wavefunction(s) for diagram number 11628
    // (none)
    // Amplitude(s) for diagram number 11628
    VVV1_0( w_fp[334], w_fp[5], w_fp[273], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];

    // *** DIAGRAM 11629 OF 15495 ***
    // Wavefunction(s) for diagram number 11629
    // (none)
    // Amplitude(s) for diagram number 11629
    FFV1_0( w_fp[3], w_fp[158], w_fp[503], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11630 OF 15495 ***
    // Wavefunction(s) for diagram number 11630
    // (none)
    // Amplitude(s) for diagram number 11630
    FFV1_0( w_fp[3], w_fp[674], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];

    // *** DIAGRAM 11631 OF 15495 ***
    // Wavefunction(s) for diagram number 11631
    // (none)
    // Amplitude(s) for diagram number 11631
    VVV1_0( w_fp[709], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    VVV1_0( w_fp[703], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    VVV1_0( w_fp[670], w_fp[164], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];

    // *** DIAGRAM 11632 OF 15495 ***
    // Wavefunction(s) for diagram number 11632
    // (none)
    // Amplitude(s) for diagram number 11632
    FFV1_0( w_fp[3], w_fp[161], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[161], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[161], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11633 OF 15495 ***
    // Wavefunction(s) for diagram number 11633
    // (none)
    // Amplitude(s) for diagram number 11633
    VVV1_0( w_fp[274], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    VVV1_0( w_fp[710], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    VVV1_0( w_fp[712], w_fp[164], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];

    // *** DIAGRAM 11634 OF 15495 ***
    // Wavefunction(s) for diagram number 11634
    // (none)
    // Amplitude(s) for diagram number 11634
    FFV1_0( w_fp[3], w_fp[158], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[158], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[158], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11635 OF 15495 ***
    // Wavefunction(s) for diagram number 11635
    // (none)
    // Amplitude(s) for diagram number 11635
    VVV1_0( w_fp[0], w_fp[341], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[342], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[343], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];

    // *** DIAGRAM 11636 OF 15495 ***
    // Wavefunction(s) for diagram number 11636
    // (none)
    // Amplitude(s) for diagram number 11636
    VVV1_0( w_fp[662], w_fp[173], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11637 OF 15495 ***
    // Wavefunction(s) for diagram number 11637
    // (none)
    // Amplitude(s) for diagram number 11637
    FFV1_0( w_fp[168], w_fp[161], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];

    // *** DIAGRAM 11638 OF 15495 ***
    // Wavefunction(s) for diagram number 11638
    // (none)
    // Amplitude(s) for diagram number 11638
    FFV1_0( w_fp[170], w_fp[156], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];

    // *** DIAGRAM 11639 OF 15495 ***
    // Wavefunction(s) for diagram number 11639
    // (none)
    // Amplitude(s) for diagram number 11639
    FFV1_0( w_fp[347], w_fp[666], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11640 OF 15495 ***
    // Wavefunction(s) for diagram number 11640
    // (none)
    // Amplitude(s) for diagram number 11640
    FFV1_0( w_fp[168], w_fp[666], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];

    // *** DIAGRAM 11641 OF 15495 ***
    // Wavefunction(s) for diagram number 11641
    // (none)
    // Amplitude(s) for diagram number 11641
    FFV1_0( w_fp[170], w_fp[666], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11642 OF 15495 ***
    // Wavefunction(s) for diagram number 11642
    // (none)
    // Amplitude(s) for diagram number 11642
    FFV1_0( w_fp[527], w_fp[506], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11643 OF 15495 ***
    // Wavefunction(s) for diagram number 11643
    // (none)
    // Amplitude(s) for diagram number 11643
    FFV1_0( w_fp[527], w_fp[156], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[286] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[346] -= amp_sv[0];

    // *** DIAGRAM 11644 OF 15495 ***
    // Wavefunction(s) for diagram number 11644
    // (none)
    // Amplitude(s) for diagram number 11644
    FFV1_0( w_fp[527], w_fp[161], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11645 OF 15495 ***
    // Wavefunction(s) for diagram number 11645
    // (none)
    // Amplitude(s) for diagram number 11645
    FFV1_0( w_fp[170], w_fp[506], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11646 OF 15495 ***
    // Wavefunction(s) for diagram number 11646
    // (none)
    // Amplitude(s) for diagram number 11646
    FFV1_0( w_fp[347], w_fp[161], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11647 OF 15495 ***
    // Wavefunction(s) for diagram number 11647
    // (none)
    // Amplitude(s) for diagram number 11647
    VVV1_0( w_fp[0], w_fp[334], w_fp[173], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11648 OF 15495 ***
    // Wavefunction(s) for diagram number 11648
    // (none)
    // Amplitude(s) for diagram number 11648
    FFV1_0( w_fp[168], w_fp[156], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[156], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[156], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11649 OF 15495 ***
    // Wavefunction(s) for diagram number 11649
    // (none)
    // Amplitude(s) for diagram number 11649
    VVV1_0( w_fp[662], w_fp[178], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11650 OF 15495 ***
    // Wavefunction(s) for diagram number 11650
    // (none)
    // Amplitude(s) for diagram number 11650
    FFV1_0( w_fp[174], w_fp[158], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];

    // *** DIAGRAM 11651 OF 15495 ***
    // Wavefunction(s) for diagram number 11651
    // (none)
    // Amplitude(s) for diagram number 11651
    FFV1_0( w_fp[175], w_fp[156], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];

    // *** DIAGRAM 11652 OF 15495 ***
    // Wavefunction(s) for diagram number 11652
    // (none)
    // Amplitude(s) for diagram number 11652
    FFV1_0( w_fp[348], w_fp[666], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11653 OF 15495 ***
    // Wavefunction(s) for diagram number 11653
    // (none)
    // Amplitude(s) for diagram number 11653
    FFV1_0( w_fp[174], w_fp[666], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];

    // *** DIAGRAM 11654 OF 15495 ***
    // Wavefunction(s) for diagram number 11654
    // (none)
    // Amplitude(s) for diagram number 11654
    FFV1_0( w_fp[175], w_fp[666], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11655 OF 15495 ***
    // Wavefunction(s) for diagram number 11655
    // (none)
    // Amplitude(s) for diagram number 11655
    FFV1_0( w_fp[594], w_fp[506], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11656 OF 15495 ***
    // Wavefunction(s) for diagram number 11656
    // (none)
    // Amplitude(s) for diagram number 11656
    FFV1_0( w_fp[594], w_fp[156], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];

    // *** DIAGRAM 11657 OF 15495 ***
    // Wavefunction(s) for diagram number 11657
    // (none)
    // Amplitude(s) for diagram number 11657
    FFV1_0( w_fp[594], w_fp[158], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11658 OF 15495 ***
    // Wavefunction(s) for diagram number 11658
    // (none)
    // Amplitude(s) for diagram number 11658
    FFV1_0( w_fp[175], w_fp[506], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11659 OF 15495 ***
    // Wavefunction(s) for diagram number 11659
    // (none)
    // Amplitude(s) for diagram number 11659
    FFV1_0( w_fp[348], w_fp[158], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11660 OF 15495 ***
    // Wavefunction(s) for diagram number 11660
    // (none)
    // Amplitude(s) for diagram number 11660
    VVV1_0( w_fp[0], w_fp[333], w_fp[178], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11661 OF 15495 ***
    // Wavefunction(s) for diagram number 11661
    // (none)
    // Amplitude(s) for diagram number 11661
    FFV1_0( w_fp[174], w_fp[156], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[156], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[156], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11662 OF 15495 ***
    // Wavefunction(s) for diagram number 11662
    // (none)
    // Amplitude(s) for diagram number 11662
    VVV1_0( w_fp[662], w_fp[164], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];

    // *** DIAGRAM 11663 OF 15495 ***
    // Wavefunction(s) for diagram number 11663
    // (none)
    // Amplitude(s) for diagram number 11663
    FFV1_0( w_fp[3], w_fp[185], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11664 OF 15495 ***
    // Wavefunction(s) for diagram number 11664
    // (none)
    // Amplitude(s) for diagram number 11664
    FFV1_0( w_fp[184], w_fp[156], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11665 OF 15495 ***
    // Wavefunction(s) for diagram number 11665
    // (none)
    // Amplitude(s) for diagram number 11665
    FFV1_0( w_fp[344], w_fp[666], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];

    // *** DIAGRAM 11666 OF 15495 ***
    // Wavefunction(s) for diagram number 11666
    // (none)
    // Amplitude(s) for diagram number 11666
    FFV1_0( w_fp[3], w_fp[666], w_fp[351], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11667 OF 15495 ***
    // Wavefunction(s) for diagram number 11667
    // (none)
    // Amplitude(s) for diagram number 11667
    FFV1_0( w_fp[184], w_fp[666], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];

    // *** DIAGRAM 11668 OF 15495 ***
    // Wavefunction(s) for diagram number 11668
    // (none)
    // Amplitude(s) for diagram number 11668
    FFV1_0( w_fp[3], w_fp[506], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11669 OF 15495 ***
    // Wavefunction(s) for diagram number 11669
    // (none)
    // Amplitude(s) for diagram number 11669
    FFV1_0( w_fp[344], w_fp[156], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11670 OF 15495 ***
    // Wavefunction(s) for diagram number 11670
    // (none)
    // Amplitude(s) for diagram number 11670
    VVV1_0( w_fp[559], w_fp[325], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];

    // *** DIAGRAM 11671 OF 15495 ***
    // Wavefunction(s) for diagram number 11671
    // (none)
    // Amplitude(s) for diagram number 11671
    FFV1_0( w_fp[184], w_fp[506], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[282] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];

    // *** DIAGRAM 11672 OF 15495 ***
    // Wavefunction(s) for diagram number 11672
    // (none)
    // Amplitude(s) for diagram number 11672
    FFV1_0( w_fp[344], w_fp[185], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];

    // *** DIAGRAM 11673 OF 15495 ***
    // Wavefunction(s) for diagram number 11673
    // (none)
    // Amplitude(s) for diagram number 11673
    VVV1_0( w_fp[0], w_fp[351], w_fp[164], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];

    // *** DIAGRAM 11674 OF 15495 ***
    // Wavefunction(s) for diagram number 11674
    // (none)
    // Amplitude(s) for diagram number 11674
    FFV1_0( w_fp[3], w_fp[156], w_fp[498], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[156], w_fp[748], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[156], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];

    // *** DIAGRAM 11675 OF 15495 ***
    // Wavefunction(s) for diagram number 11675
    // (none)
    // Amplitude(s) for diagram number 11675
    VVVV1_0( w_fp[662], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[194], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 11676 OF 15495 ***
    // Wavefunction(s) for diagram number 11676
    // (none)
    // Amplitude(s) for diagram number 11676
    VVV1_0( w_fp[194], w_fp[6], w_fp[548], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];

    // *** DIAGRAM 11677 OF 15495 ***
    // Wavefunction(s) for diagram number 11677
    // (none)
    // Amplitude(s) for diagram number 11677
    VVV1_0( w_fp[194], w_fp[4], w_fp[500], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 11678 OF 15495 ***
    // Wavefunction(s) for diagram number 11678
    // (none)
    // Amplitude(s) for diagram number 11678
    FFV1_0( w_fp[554], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];

    // *** DIAGRAM 11679 OF 15495 ***
    // Wavefunction(s) for diagram number 11679
    // (none)
    // Amplitude(s) for diagram number 11679
    FFV1_0( w_fp[3], w_fp[190], w_fp[500], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11680 OF 15495 ***
    // Wavefunction(s) for diagram number 11680
    // (none)
    // Amplitude(s) for diagram number 11680
    FFV1_0( w_fp[554], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];

    // *** DIAGRAM 11681 OF 15495 ***
    // Wavefunction(s) for diagram number 11681
    // (none)
    // Amplitude(s) for diagram number 11681
    FFV1_0( w_fp[3], w_fp[191], w_fp[548], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11682 OF 15495 ***
    // Wavefunction(s) for diagram number 11682
    // (none)
    // Amplitude(s) for diagram number 11682
    FFV1_0( w_fp[344], w_fp[679], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11683 OF 15495 ***
    // Wavefunction(s) for diagram number 11683
    // (none)
    // Amplitude(s) for diagram number 11683
    FFV1_0( w_fp[344], w_fp[677], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11684 OF 15495 ***
    // Wavefunction(s) for diagram number 11684
    // (none)
    // Amplitude(s) for diagram number 11684
    VVV1_0( w_fp[332], w_fp[6], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11685 OF 15495 ***
    // Wavefunction(s) for diagram number 11685
    // (none)
    // Amplitude(s) for diagram number 11685
    FFV1_0( w_fp[3], w_fp[677], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];

    // *** DIAGRAM 11686 OF 15495 ***
    // Wavefunction(s) for diagram number 11686
    // (none)
    // Amplitude(s) for diagram number 11686
    VVV1_0( w_fp[334], w_fp[4], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11687 OF 15495 ***
    // Wavefunction(s) for diagram number 11687
    // (none)
    // Amplitude(s) for diagram number 11687
    FFV1_0( w_fp[3], w_fp[679], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[367] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];

    // *** DIAGRAM 11688 OF 15495 ***
    // Wavefunction(s) for diagram number 11688
    // (none)
    // Amplitude(s) for diagram number 11688
    FFV1_0( w_fp[3], w_fp[676], w_fp[338], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[339], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[340], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11689 OF 15495 ***
    // Wavefunction(s) for diagram number 11689
    // (none)
    // Amplitude(s) for diagram number 11689
    FFV1_0( w_fp[317], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11690 OF 15495 ***
    // Wavefunction(s) for diagram number 11690
    // (none)
    // Amplitude(s) for diagram number 11690
    FFV1_0( w_fp[344], w_fp[683], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11691 OF 15495 ***
    // Wavefunction(s) for diagram number 11691
    // (none)
    // Amplitude(s) for diagram number 11691
    FFV1_0( w_fp[317], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11692 OF 15495 ***
    // Wavefunction(s) for diagram number 11692
    // (none)
    // Amplitude(s) for diagram number 11692
    FFV1_0( w_fp[344], w_fp[681], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11693 OF 15495 ***
    // Wavefunction(s) for diagram number 11693
    // (none)
    // Amplitude(s) for diagram number 11693
    VVVV1_0( w_fp[0], w_fp[332], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[332], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[332], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];

    // *** DIAGRAM 11694 OF 15495 ***
    // Wavefunction(s) for diagram number 11694
    // (none)
    // Amplitude(s) for diagram number 11694
    VVV1_0( w_fp[194], w_fp[6], w_fp[743], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];

    // *** DIAGRAM 11695 OF 15495 ***
    // Wavefunction(s) for diagram number 11695
    // (none)
    // Amplitude(s) for diagram number 11695
    VVV1_0( w_fp[332], w_fp[6], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];

    // *** DIAGRAM 11696 OF 15495 ***
    // Wavefunction(s) for diagram number 11696
    // (none)
    // Amplitude(s) for diagram number 11696
    FFV1_0( w_fp[3], w_fp[191], w_fp[743], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11697 OF 15495 ***
    // Wavefunction(s) for diagram number 11697
    // (none)
    // Amplitude(s) for diagram number 11697
    FFV1_0( w_fp[3], w_fp[681], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];

    // *** DIAGRAM 11698 OF 15495 ***
    // Wavefunction(s) for diagram number 11698
    // (none)
    // Amplitude(s) for diagram number 11698
    VVVV1_0( w_fp[0], w_fp[334], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[334], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[334], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];

    // *** DIAGRAM 11699 OF 15495 ***
    // Wavefunction(s) for diagram number 11699
    // (none)
    // Amplitude(s) for diagram number 11699
    VVV1_0( w_fp[194], w_fp[4], w_fp[503], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];

    // *** DIAGRAM 11700 OF 15495 ***
    // Wavefunction(s) for diagram number 11700
    // (none)
    // Amplitude(s) for diagram number 11700
    VVV1_0( w_fp[334], w_fp[4], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];

    // *** DIAGRAM 11701 OF 15495 ***
    // Wavefunction(s) for diagram number 11701
    // (none)
    // Amplitude(s) for diagram number 11701
    FFV1_0( w_fp[3], w_fp[190], w_fp[503], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11702 OF 15495 ***
    // Wavefunction(s) for diagram number 11702
    // (none)
    // Amplitude(s) for diagram number 11702
    FFV1_0( w_fp[3], w_fp[683], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];

    // *** DIAGRAM 11703 OF 15495 ***
    // Wavefunction(s) for diagram number 11703
    // (none)
    // Amplitude(s) for diagram number 11703
    VVV1_0( w_fp[551], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    VVV1_0( w_fp[504], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    VVV1_0( w_fp[544], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];

    // *** DIAGRAM 11704 OF 15495 ***
    // Wavefunction(s) for diagram number 11704
    // (none)
    // Amplitude(s) for diagram number 11704
    FFV1_0( w_fp[3], w_fp[191], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[191], w_fp[504], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[191], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11705 OF 15495 ***
    // Wavefunction(s) for diagram number 11705
    // (none)
    // Amplitude(s) for diagram number 11705
    VVV1_0( w_fp[274], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    VVV1_0( w_fp[710], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    VVV1_0( w_fp[712], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];

    // *** DIAGRAM 11706 OF 15495 ***
    // Wavefunction(s) for diagram number 11706
    // (none)
    // Amplitude(s) for diagram number 11706
    FFV1_0( w_fp[3], w_fp[190], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[190], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[190], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11707 OF 15495 ***
    // Wavefunction(s) for diagram number 11707
    // (none)
    // Amplitude(s) for diagram number 11707
    VVV1_0( w_fp[0], w_fp[338], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[339], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[340], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];

    // *** DIAGRAM 11708 OF 15495 ***
    // Wavefunction(s) for diagram number 11708
    // (none)
    // Amplitude(s) for diagram number 11708
    VVV1_0( w_fp[662], w_fp[201], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11709 OF 15495 ***
    // Wavefunction(s) for diagram number 11709
    // (none)
    // Amplitude(s) for diagram number 11709
    FFV1_0( w_fp[196], w_fp[191], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];

    // *** DIAGRAM 11710 OF 15495 ***
    // Wavefunction(s) for diagram number 11710
    // (none)
    // Amplitude(s) for diagram number 11710
    FFV1_0( w_fp[198], w_fp[169], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 11711 OF 15495 ***
    // Wavefunction(s) for diagram number 11711
    // (none)
    // Amplitude(s) for diagram number 11711
    FFV1_0( w_fp[346], w_fp[676], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11712 OF 15495 ***
    // Wavefunction(s) for diagram number 11712
    // (none)
    // Amplitude(s) for diagram number 11712
    FFV1_0( w_fp[196], w_fp[676], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];

    // *** DIAGRAM 11713 OF 15495 ***
    // Wavefunction(s) for diagram number 11713
    // (none)
    // Amplitude(s) for diagram number 11713
    FFV1_0( w_fp[198], w_fp[676], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11714 OF 15495 ***
    // Wavefunction(s) for diagram number 11714
    // (none)
    // Amplitude(s) for diagram number 11714
    FFV1_0( w_fp[536], w_fp[507], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11715 OF 15495 ***
    // Wavefunction(s) for diagram number 11715
    // (none)
    // Amplitude(s) for diagram number 11715
    FFV1_0( w_fp[536], w_fp[169], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[406] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];

    // *** DIAGRAM 11716 OF 15495 ***
    // Wavefunction(s) for diagram number 11716
    // (none)
    // Amplitude(s) for diagram number 11716
    FFV1_0( w_fp[536], w_fp[191], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11717 OF 15495 ***
    // Wavefunction(s) for diagram number 11717
    // (none)
    // Amplitude(s) for diagram number 11717
    FFV1_0( w_fp[198], w_fp[507], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11718 OF 15495 ***
    // Wavefunction(s) for diagram number 11718
    // (none)
    // Amplitude(s) for diagram number 11718
    FFV1_0( w_fp[346], w_fp[191], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11719 OF 15495 ***
    // Wavefunction(s) for diagram number 11719
    // (none)
    // Amplitude(s) for diagram number 11719
    VVV1_0( w_fp[0], w_fp[334], w_fp[201], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11720 OF 15495 ***
    // Wavefunction(s) for diagram number 11720
    // (none)
    // Amplitude(s) for diagram number 11720
    FFV1_0( w_fp[196], w_fp[169], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[169], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[169], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11721 OF 15495 ***
    // Wavefunction(s) for diagram number 11721
    // (none)
    // Amplitude(s) for diagram number 11721
    VVV1_0( w_fp[662], w_fp[203], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11722 OF 15495 ***
    // Wavefunction(s) for diagram number 11722
    // (none)
    // Amplitude(s) for diagram number 11722
    FFV1_0( w_fp[174], w_fp[190], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];

    // *** DIAGRAM 11723 OF 15495 ***
    // Wavefunction(s) for diagram number 11723
    // (none)
    // Amplitude(s) for diagram number 11723
    FFV1_0( w_fp[202], w_fp[169], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];

    // *** DIAGRAM 11724 OF 15495 ***
    // Wavefunction(s) for diagram number 11724
    // (none)
    // Amplitude(s) for diagram number 11724
    FFV1_0( w_fp[348], w_fp[676], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11725 OF 15495 ***
    // Wavefunction(s) for diagram number 11725
    // (none)
    // Amplitude(s) for diagram number 11725
    FFV1_0( w_fp[174], w_fp[676], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];

    // *** DIAGRAM 11726 OF 15495 ***
    // Wavefunction(s) for diagram number 11726
    // (none)
    // Amplitude(s) for diagram number 11726
    FFV1_0( w_fp[202], w_fp[676], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11727 OF 15495 ***
    // Wavefunction(s) for diagram number 11727
    // (none)
    // Amplitude(s) for diagram number 11727
    FFV1_0( w_fp[594], w_fp[507], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11728 OF 15495 ***
    // Wavefunction(s) for diagram number 11728
    // (none)
    // Amplitude(s) for diagram number 11728
    FFV1_0( w_fp[594], w_fp[169], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];

    // *** DIAGRAM 11729 OF 15495 ***
    // Wavefunction(s) for diagram number 11729
    // (none)
    // Amplitude(s) for diagram number 11729
    FFV1_0( w_fp[594], w_fp[190], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11730 OF 15495 ***
    // Wavefunction(s) for diagram number 11730
    // (none)
    // Amplitude(s) for diagram number 11730
    FFV1_0( w_fp[202], w_fp[507], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11731 OF 15495 ***
    // Wavefunction(s) for diagram number 11731
    // (none)
    // Amplitude(s) for diagram number 11731
    FFV1_0( w_fp[348], w_fp[190], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11732 OF 15495 ***
    // Wavefunction(s) for diagram number 11732
    // (none)
    // Amplitude(s) for diagram number 11732
    VVV1_0( w_fp[0], w_fp[332], w_fp[203], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11733 OF 15495 ***
    // Wavefunction(s) for diagram number 11733
    // (none)
    // Amplitude(s) for diagram number 11733
    FFV1_0( w_fp[174], w_fp[169], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[504], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11734 OF 15495 ***
    // Wavefunction(s) for diagram number 11734
    // (none)
    // Amplitude(s) for diagram number 11734
    VVV1_0( w_fp[662], w_fp[194], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];

    // *** DIAGRAM 11735 OF 15495 ***
    // Wavefunction(s) for diagram number 11735
    // (none)
    // Amplitude(s) for diagram number 11735
    FFV1_0( w_fp[3], w_fp[207], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11736 OF 15495 ***
    // Wavefunction(s) for diagram number 11736
    // (none)
    // Amplitude(s) for diagram number 11736
    FFV1_0( w_fp[206], w_fp[169], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11737 OF 15495 ***
    // Wavefunction(s) for diagram number 11737
    // (none)
    // Amplitude(s) for diagram number 11737
    FFV1_0( w_fp[344], w_fp[676], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];

    // *** DIAGRAM 11738 OF 15495 ***
    // Wavefunction(s) for diagram number 11738
    // (none)
    // Amplitude(s) for diagram number 11738
    FFV1_0( w_fp[3], w_fp[676], w_fp[350], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11739 OF 15495 ***
    // Wavefunction(s) for diagram number 11739
    // (none)
    // Amplitude(s) for diagram number 11739
    FFV1_0( w_fp[206], w_fp[676], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];

    // *** DIAGRAM 11740 OF 15495 ***
    // Wavefunction(s) for diagram number 11740
    // (none)
    // Amplitude(s) for diagram number 11740
    FFV1_0( w_fp[3], w_fp[507], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11741 OF 15495 ***
    // Wavefunction(s) for diagram number 11741
    // (none)
    // Amplitude(s) for diagram number 11741
    FFV1_0( w_fp[344], w_fp[169], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11742 OF 15495 ***
    // Wavefunction(s) for diagram number 11742
    // (none)
    // Amplitude(s) for diagram number 11742
    VVV1_0( w_fp[576], w_fp[325], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];

    // *** DIAGRAM 11743 OF 15495 ***
    // Wavefunction(s) for diagram number 11743
    // (none)
    // Amplitude(s) for diagram number 11743
    FFV1_0( w_fp[206], w_fp[507], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 11744 OF 15495 ***
    // Wavefunction(s) for diagram number 11744
    // (none)
    // Amplitude(s) for diagram number 11744
    FFV1_0( w_fp[344], w_fp[207], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];

    // *** DIAGRAM 11745 OF 15495 ***
    // Wavefunction(s) for diagram number 11745
    // (none)
    // Amplitude(s) for diagram number 11745
    VVV1_0( w_fp[0], w_fp[350], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];

    // *** DIAGRAM 11746 OF 15495 ***
    // Wavefunction(s) for diagram number 11746
    // (none)
    // Amplitude(s) for diagram number 11746
    FFV1_0( w_fp[3], w_fp[169], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[329], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[328], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];

    // *** DIAGRAM 11747 OF 15495 ***
    // Wavefunction(s) for diagram number 11747
    // (none)
    // Amplitude(s) for diagram number 11747
    VVVV1_0( w_fp[662], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 11748 OF 15495 ***
    // Wavefunction(s) for diagram number 11748
    // (none)
    // Amplitude(s) for diagram number 11748
    VVV1_0( w_fp[214], w_fp[5], w_fp[548], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];

    // *** DIAGRAM 11749 OF 15495 ***
    // Wavefunction(s) for diagram number 11749
    // (none)
    // Amplitude(s) for diagram number 11749
    VVV1_0( w_fp[214], w_fp[4], w_fp[706], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 11750 OF 15495 ***
    // Wavefunction(s) for diagram number 11750
    // (none)
    // Amplitude(s) for diagram number 11750
    FFV1_0( w_fp[554], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];

    // *** DIAGRAM 11751 OF 15495 ***
    // Wavefunction(s) for diagram number 11751
    // (none)
    // Amplitude(s) for diagram number 11751
    FFV1_0( w_fp[3], w_fp[211], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11752 OF 15495 ***
    // Wavefunction(s) for diagram number 11752
    // (none)
    // Amplitude(s) for diagram number 11752
    FFV1_0( w_fp[554], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 11753 OF 15495 ***
    // Wavefunction(s) for diagram number 11753
    // (none)
    // Amplitude(s) for diagram number 11753
    FFV1_0( w_fp[3], w_fp[212], w_fp[548], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11754 OF 15495 ***
    // Wavefunction(s) for diagram number 11754
    // (none)
    // Amplitude(s) for diagram number 11754
    FFV1_0( w_fp[344], w_fp[686], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11755 OF 15495 ***
    // Wavefunction(s) for diagram number 11755
    // (none)
    // Amplitude(s) for diagram number 11755
    FFV1_0( w_fp[344], w_fp[684], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11756 OF 15495 ***
    // Wavefunction(s) for diagram number 11756
    // (none)
    // Amplitude(s) for diagram number 11756
    VVV1_0( w_fp[332], w_fp[5], w_fp[489], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11757 OF 15495 ***
    // Wavefunction(s) for diagram number 11757
    // (none)
    // Amplitude(s) for diagram number 11757
    FFV1_0( w_fp[3], w_fp[684], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[493] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];

    // *** DIAGRAM 11758 OF 15495 ***
    // Wavefunction(s) for diagram number 11758
    // (none)
    // Amplitude(s) for diagram number 11758
    VVV1_0( w_fp[333], w_fp[4], w_fp[489], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11759 OF 15495 ***
    // Wavefunction(s) for diagram number 11759
    // (none)
    // Amplitude(s) for diagram number 11759
    FFV1_0( w_fp[3], w_fp[686], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[487] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];

    // *** DIAGRAM 11760 OF 15495 ***
    // Wavefunction(s) for diagram number 11760
    // (none)
    // Amplitude(s) for diagram number 11760
    FFV1_0( w_fp[3], w_fp[540], w_fp[335], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[540], w_fp[336], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[540], w_fp[337], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11761 OF 15495 ***
    // Wavefunction(s) for diagram number 11761
    // (none)
    // Amplitude(s) for diagram number 11761
    FFV1_0( w_fp[317], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11762 OF 15495 ***
    // Wavefunction(s) for diagram number 11762
    // (none)
    // Amplitude(s) for diagram number 11762
    FFV1_0( w_fp[344], w_fp[689], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11763 OF 15495 ***
    // Wavefunction(s) for diagram number 11763
    // (none)
    // Amplitude(s) for diagram number 11763
    FFV1_0( w_fp[317], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11764 OF 15495 ***
    // Wavefunction(s) for diagram number 11764
    // (none)
    // Amplitude(s) for diagram number 11764
    FFV1_0( w_fp[344], w_fp[687], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11765 OF 15495 ***
    // Wavefunction(s) for diagram number 11765
    // (none)
    // Amplitude(s) for diagram number 11765
    VVVV1_0( w_fp[0], w_fp[332], w_fp[214], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[332], w_fp[214], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[332], w_fp[214], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];

    // *** DIAGRAM 11766 OF 15495 ***
    // Wavefunction(s) for diagram number 11766
    // (none)
    // Amplitude(s) for diagram number 11766
    VVV1_0( w_fp[214], w_fp[5], w_fp[743], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];

    // *** DIAGRAM 11767 OF 15495 ***
    // Wavefunction(s) for diagram number 11767
    // (none)
    // Amplitude(s) for diagram number 11767
    VVV1_0( w_fp[332], w_fp[5], w_fp[745], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];

    // *** DIAGRAM 11768 OF 15495 ***
    // Wavefunction(s) for diagram number 11768
    // (none)
    // Amplitude(s) for diagram number 11768
    FFV1_0( w_fp[3], w_fp[212], w_fp[743], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11769 OF 15495 ***
    // Wavefunction(s) for diagram number 11769
    // (none)
    // Amplitude(s) for diagram number 11769
    FFV1_0( w_fp[3], w_fp[687], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];

    // *** DIAGRAM 11770 OF 15495 ***
    // Wavefunction(s) for diagram number 11770
    // (none)
    // Amplitude(s) for diagram number 11770
    VVVV1_0( w_fp[0], w_fp[333], w_fp[214], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[333], w_fp[214], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[333], w_fp[214], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];

    // *** DIAGRAM 11771 OF 15495 ***
    // Wavefunction(s) for diagram number 11771
    // (none)
    // Amplitude(s) for diagram number 11771
    VVV1_0( w_fp[214], w_fp[4], w_fp[270], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];

    // *** DIAGRAM 11772 OF 15495 ***
    // Wavefunction(s) for diagram number 11772
    // (none)
    // Amplitude(s) for diagram number 11772
    VVV1_0( w_fp[333], w_fp[4], w_fp[745], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];

    // *** DIAGRAM 11773 OF 15495 ***
    // Wavefunction(s) for diagram number 11773
    // (none)
    // Amplitude(s) for diagram number 11773
    FFV1_0( w_fp[3], w_fp[211], w_fp[270], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11774 OF 15495 ***
    // Wavefunction(s) for diagram number 11774
    // (none)
    // Amplitude(s) for diagram number 11774
    FFV1_0( w_fp[3], w_fp[689], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];

    // *** DIAGRAM 11775 OF 15495 ***
    // Wavefunction(s) for diagram number 11775
    // (none)
    // Amplitude(s) for diagram number 11775
    VVV1_0( w_fp[551], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    VVV1_0( w_fp[504], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    VVV1_0( w_fp[544], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];

    // *** DIAGRAM 11776 OF 15495 ***
    // Wavefunction(s) for diagram number 11776
    // (none)
    // Amplitude(s) for diagram number 11776
    FFV1_0( w_fp[3], w_fp[212], w_fp[551], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[212], w_fp[504], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[212], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11777 OF 15495 ***
    // Wavefunction(s) for diagram number 11777
    // (none)
    // Amplitude(s) for diagram number 11777
    VVV1_0( w_fp[709], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    VVV1_0( w_fp[703], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    VVV1_0( w_fp[670], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];

    // *** DIAGRAM 11778 OF 15495 ***
    // Wavefunction(s) for diagram number 11778
    // (none)
    // Amplitude(s) for diagram number 11778
    FFV1_0( w_fp[3], w_fp[211], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[211], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[211], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11779 OF 15495 ***
    // Wavefunction(s) for diagram number 11779
    // (none)
    // Amplitude(s) for diagram number 11779
    VVV1_0( w_fp[0], w_fp[335], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[336], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[337], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];

    // *** DIAGRAM 11780 OF 15495 ***
    // Wavefunction(s) for diagram number 11780
    // (none)
    // Amplitude(s) for diagram number 11780
    VVV1_0( w_fp[662], w_fp[217], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11781 OF 15495 ***
    // Wavefunction(s) for diagram number 11781
    // (none)
    // Amplitude(s) for diagram number 11781
    FFV1_0( w_fp[196], w_fp[212], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];

    // *** DIAGRAM 11782 OF 15495 ***
    // Wavefunction(s) for diagram number 11782
    // (none)
    // Amplitude(s) for diagram number 11782
    FFV1_0( w_fp[216], w_fp[197], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 11783 OF 15495 ***
    // Wavefunction(s) for diagram number 11783
    // (none)
    // Amplitude(s) for diagram number 11783
    FFV1_0( w_fp[346], w_fp[540], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11784 OF 15495 ***
    // Wavefunction(s) for diagram number 11784
    // (none)
    // Amplitude(s) for diagram number 11784
    FFV1_0( w_fp[196], w_fp[540], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];

    // *** DIAGRAM 11785 OF 15495 ***
    // Wavefunction(s) for diagram number 11785
    // (none)
    // Amplitude(s) for diagram number 11785
    FFV1_0( w_fp[216], w_fp[540], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11786 OF 15495 ***
    // Wavefunction(s) for diagram number 11786
    // (none)
    // Amplitude(s) for diagram number 11786
    FFV1_0( w_fp[536], w_fp[508], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11787 OF 15495 ***
    // Wavefunction(s) for diagram number 11787
    // (none)
    // Amplitude(s) for diagram number 11787
    FFV1_0( w_fp[536], w_fp[197], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[526] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];

    // *** DIAGRAM 11788 OF 15495 ***
    // Wavefunction(s) for diagram number 11788
    // (none)
    // Amplitude(s) for diagram number 11788
    FFV1_0( w_fp[536], w_fp[212], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11789 OF 15495 ***
    // Wavefunction(s) for diagram number 11789
    // (none)
    // Amplitude(s) for diagram number 11789
    FFV1_0( w_fp[216], w_fp[508], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11790 OF 15495 ***
    // Wavefunction(s) for diagram number 11790
    // (none)
    // Amplitude(s) for diagram number 11790
    FFV1_0( w_fp[346], w_fp[212], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11791 OF 15495 ***
    // Wavefunction(s) for diagram number 11791
    // (none)
    // Amplitude(s) for diagram number 11791
    VVV1_0( w_fp[0], w_fp[333], w_fp[217], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11792 OF 15495 ***
    // Wavefunction(s) for diagram number 11792
    // (none)
    // Amplitude(s) for diagram number 11792
    FFV1_0( w_fp[196], w_fp[197], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[197], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[197], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11793 OF 15495 ***
    // Wavefunction(s) for diagram number 11793
    // (none)
    // Amplitude(s) for diagram number 11793
    VVV1_0( w_fp[662], w_fp[219], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11794 OF 15495 ***
    // Wavefunction(s) for diagram number 11794
    // (none)
    // Amplitude(s) for diagram number 11794
    FFV1_0( w_fp[168], w_fp[211], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];

    // *** DIAGRAM 11795 OF 15495 ***
    // Wavefunction(s) for diagram number 11795
    // (none)
    // Amplitude(s) for diagram number 11795
    FFV1_0( w_fp[218], w_fp[197], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];

    // *** DIAGRAM 11796 OF 15495 ***
    // Wavefunction(s) for diagram number 11796
    // (none)
    // Amplitude(s) for diagram number 11796
    FFV1_0( w_fp[347], w_fp[540], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11797 OF 15495 ***
    // Wavefunction(s) for diagram number 11797
    // (none)
    // Amplitude(s) for diagram number 11797
    FFV1_0( w_fp[168], w_fp[540], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];

    // *** DIAGRAM 11798 OF 15495 ***
    // Wavefunction(s) for diagram number 11798
    // (none)
    // Amplitude(s) for diagram number 11798
    FFV1_0( w_fp[218], w_fp[540], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11799 OF 15495 ***
    // Wavefunction(s) for diagram number 11799
    // (none)
    // Amplitude(s) for diagram number 11799
    FFV1_0( w_fp[527], w_fp[508], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11800 OF 15495 ***
    // Wavefunction(s) for diagram number 11800
    // (none)
    // Amplitude(s) for diagram number 11800
    FFV1_0( w_fp[527], w_fp[197], w_fp[332], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[524] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 317 );
    storeWf( wfs, w_cx, nevt, 554 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
