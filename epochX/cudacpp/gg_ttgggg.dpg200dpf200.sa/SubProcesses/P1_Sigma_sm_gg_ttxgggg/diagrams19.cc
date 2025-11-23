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
  diagramgroup19( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 93 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 119 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 142 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 165 );
    retrieveWf( wfs, w_cx, nevt, 166 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 172 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 195 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 229 );
    retrieveWf( wfs, w_cx, nevt, 240 );
    retrieveWf( wfs, w_cx, nevt, 241 );
    retrieveWf( wfs, w_cx, nevt, 242 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 276 );
    retrieveWf( wfs, w_cx, nevt, 277 );
    retrieveWf( wfs, w_cx, nevt, 280 );
    retrieveWf( wfs, w_cx, nevt, 282 );
    retrieveWf( wfs, w_cx, nevt, 284 );
    retrieveWf( wfs, w_cx, nevt, 303 );
    retrieveWf( wfs, w_cx, nevt, 305 );
    retrieveWf( wfs, w_cx, nevt, 345 );
    retrieveWf( wfs, w_cx, nevt, 356 );
    retrieveWf( wfs, w_cx, nevt, 357 );
    retrieveWf( wfs, w_cx, nevt, 378 );
    retrieveWf( wfs, w_cx, nevt, 379 );
    retrieveWf( wfs, w_cx, nevt, 380 );
    retrieveWf( wfs, w_cx, nevt, 381 );
    retrieveWf( wfs, w_cx, nevt, 382 );
    retrieveWf( wfs, w_cx, nevt, 383 );
    retrieveWf( wfs, w_cx, nevt, 384 );
    retrieveWf( wfs, w_cx, nevt, 385 );
    retrieveWf( wfs, w_cx, nevt, 389 );
    retrieveWf( wfs, w_cx, nevt, 390 );
    retrieveWf( wfs, w_cx, nevt, 391 );
    retrieveWf( wfs, w_cx, nevt, 392 );
    retrieveWf( wfs, w_cx, nevt, 393 );
    retrieveWf( wfs, w_cx, nevt, 394 );
    retrieveWf( wfs, w_cx, nevt, 395 );
    retrieveWf( wfs, w_cx, nevt, 396 );
    retrieveWf( wfs, w_cx, nevt, 400 );
    retrieveWf( wfs, w_cx, nevt, 401 );
    retrieveWf( wfs, w_cx, nevt, 402 );
    retrieveWf( wfs, w_cx, nevt, 403 );
    retrieveWf( wfs, w_cx, nevt, 404 );
    retrieveWf( wfs, w_cx, nevt, 405 );
    retrieveWf( wfs, w_cx, nevt, 406 );
    retrieveWf( wfs, w_cx, nevt, 407 );
    retrieveWf( wfs, w_cx, nevt, 411 );
    retrieveWf( wfs, w_cx, nevt, 412 );
    retrieveWf( wfs, w_cx, nevt, 413 );
    retrieveWf( wfs, w_cx, nevt, 414 );
    retrieveWf( wfs, w_cx, nevt, 415 );
    retrieveWf( wfs, w_cx, nevt, 416 );
    retrieveWf( wfs, w_cx, nevt, 417 );
    retrieveWf( wfs, w_cx, nevt, 421 );
    retrieveWf( wfs, w_cx, nevt, 422 );
    retrieveWf( wfs, w_cx, nevt, 423 );
    retrieveWf( wfs, w_cx, nevt, 424 );
    retrieveWf( wfs, w_cx, nevt, 425 );
    retrieveWf( wfs, w_cx, nevt, 426 );
    retrieveWf( wfs, w_cx, nevt, 427 );
    retrieveWf( wfs, w_cx, nevt, 431 );
    retrieveWf( wfs, w_cx, nevt, 432 );
    retrieveWf( wfs, w_cx, nevt, 433 );
    retrieveWf( wfs, w_cx, nevt, 434 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 439 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 446 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 459 );
    retrieveWf( wfs, w_cx, nevt, 461 );
    retrieveWf( wfs, w_cx, nevt, 462 );
    retrieveWf( wfs, w_cx, nevt, 464 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 484 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 492 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 513 );
#endif
#endif

    // *** DIAGRAM 3601 OF 15495 ***
    // Wavefunction(s) for diagram number 3601
    // (none)
    // Amplitude(s) for diagram number 3601
    FFV1_0( w_fp[447], w_fp[2], w_fp[277], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[357], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[74], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3602 OF 15495 ***
    // Wavefunction(s) for diagram number 3602
    // (none)
    // Amplitude(s) for diagram number 3602
    FFV1_0( w_fp[451], w_fp[484], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[513], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[446], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];

    // *** DIAGRAM 3603 OF 15495 ***
    // Wavefunction(s) for diagram number 3603
    // (none)
    // Amplitude(s) for diagram number 3603
    FFV1_0( w_fp[451], w_fp[2], w_fp[303], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[172], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[76], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3604 OF 15495 ***
    // Wavefunction(s) for diagram number 3604
    FFV1_2( w_fp[45], w_fp[276], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    FFV1_2( w_fp[45], w_fp[356], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[455] );
    FFV1_2( w_fp[45], w_fp[142], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[440] );
    // Amplitude(s) for diagram number 3604
    FFV1_0( w_fp[495], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    FFV1_0( w_fp[455], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    FFV1_0( w_fp[440], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 3605 OF 15495 ***
    // Wavefunction(s) for diagram number 3605
    // (none)
    // Amplitude(s) for diagram number 3605
    VVV1_0( w_fp[276], w_fp[7], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[356], w_fp[7], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[142], w_fp[7], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3606 OF 15495 ***
    // Wavefunction(s) for diagram number 3606
    // (none)
    // Amplitude(s) for diagram number 3606
    FFV1_0( w_fp[451], w_fp[197], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[197], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[197], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];

    // *** DIAGRAM 3607 OF 15495 ***
    // Wavefunction(s) for diagram number 3607
    // (none)
    // Amplitude(s) for diagram number 3607
    FFV1_0( w_fp[495], w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    FFV1_0( w_fp[455], w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    FFV1_0( w_fp[440], w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3608 OF 15495 ***
    // Wavefunction(s) for diagram number 3608
    // (none)
    // Amplitude(s) for diagram number 3608
    VVV1_0( w_fp[276], w_fp[6], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[356], w_fp[6], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[142], w_fp[6], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3609 OF 15495 ***
    // Wavefunction(s) for diagram number 3609
    // (none)
    // Amplitude(s) for diagram number 3609
    FFV1_0( w_fp[447], w_fp[215], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[215], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[215], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 3610 OF 15495 ***
    // Wavefunction(s) for diagram number 3610
    // (none)
    // Amplitude(s) for diagram number 3610
    FFV1_0( w_fp[495], w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[455], w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[440], w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3611 OF 15495 ***
    // Wavefunction(s) for diagram number 3611
    // (none)
    // Amplitude(s) for diagram number 3611
    VVV1_0( w_fp[276], w_fp[84], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[356], w_fp[84], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVV1_0( w_fp[142], w_fp[84], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3612 OF 15495 ***
    // Wavefunction(s) for diagram number 3612
    // (none)
    // Amplitude(s) for diagram number 3612
    FFV1_0( w_fp[458], w_fp[2], w_fp[276], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[458], w_fp[2], w_fp[356], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[458], w_fp[2], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3613 OF 15495 ***
    // Wavefunction(s) for diagram number 3613
    // (none)
    // Amplitude(s) for diagram number 3613
    VVV1_0( w_fp[479], w_fp[378], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[381], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[382], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3614 OF 15495 ***
    // Wavefunction(s) for diagram number 3614
    // (none)
    // Amplitude(s) for diagram number 3614
    VVV1_0( w_fp[479], w_fp[383], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[384], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[385], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3615 OF 15495 ***
    // Wavefunction(s) for diagram number 3615
    // (none)
    // Amplitude(s) for diagram number 3615
    VVVV1_0( w_fp[93], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVVV3_0( w_fp[93], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVVV4_0( w_fp[93], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    VVVV1_0( w_fp[380], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    VVVV3_0( w_fp[380], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    VVVV4_0( w_fp[380], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    VVVV1_0( w_fp[379], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVVV3_0( w_fp[379], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVVV4_0( w_fp[379], w_fp[5], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 3616 OF 15495 ***
    // Wavefunction(s) for diagram number 3616
    FFV1_1( w_fp[2], w_fp[93], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[458] );
    FFV1_1( w_fp[2], w_fp[380], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[440] );
    FFV1_1( w_fp[2], w_fp[379], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[455] );
    // Amplitude(s) for diagram number 3616
    FFV1_0( w_fp[442], w_fp[458], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[440], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[455], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];

    // *** DIAGRAM 3617 OF 15495 ***
    // Wavefunction(s) for diagram number 3617
    // (none)
    // Amplitude(s) for diagram number 3617
    FFV1_0( w_fp[442], w_fp[2], w_fp[383], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[384], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[385], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3618 OF 15495 ***
    // Wavefunction(s) for diagram number 3618
    // (none)
    // Amplitude(s) for diagram number 3618
    FFV1_0( w_fp[451], w_fp[458], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[440], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[455], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];

    // *** DIAGRAM 3619 OF 15495 ***
    // Wavefunction(s) for diagram number 3619
    // (none)
    // Amplitude(s) for diagram number 3619
    FFV1_0( w_fp[451], w_fp[2], w_fp[378], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[381], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[382], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3620 OF 15495 ***
    // Wavefunction(s) for diagram number 3620
    FFV1_2( w_fp[45], w_fp[93], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    FFV1_2( w_fp[45], w_fp[380], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[487] );
    FFV1_2( w_fp[45], w_fp[379], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[469] );
    // Amplitude(s) for diagram number 3620
    FFV1_0( w_fp[495], w_fp[169], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    FFV1_0( w_fp[487], w_fp[169], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    FFV1_0( w_fp[469], w_fp[169], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 3621 OF 15495 ***
    // Wavefunction(s) for diagram number 3621
    // (none)
    // Amplitude(s) for diagram number 3621
    VVV1_0( w_fp[93], w_fp[7], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[380], w_fp[7], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[379], w_fp[7], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3622 OF 15495 ***
    // Wavefunction(s) for diagram number 3622
    // (none)
    // Amplitude(s) for diagram number 3622
    FFV1_0( w_fp[451], w_fp[169], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[169], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[169], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];

    // *** DIAGRAM 3623 OF 15495 ***
    // Wavefunction(s) for diagram number 3623
    // (none)
    // Amplitude(s) for diagram number 3623
    FFV1_0( w_fp[495], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    FFV1_0( w_fp[487], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    FFV1_0( w_fp[469], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 3624 OF 15495 ***
    // Wavefunction(s) for diagram number 3624
    // (none)
    // Amplitude(s) for diagram number 3624
    VVV1_0( w_fp[93], w_fp[5], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[380], w_fp[5], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[379], w_fp[5], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3625 OF 15495 ***
    // Wavefunction(s) for diagram number 3625
    // (none)
    // Amplitude(s) for diagram number 3625
    FFV1_0( w_fp[442], w_fp[215], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[215], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[215], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3626 OF 15495 ***
    // Wavefunction(s) for diagram number 3626
    // (none)
    // Amplitude(s) for diagram number 3626
    FFV1_0( w_fp[495], w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[487], w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[469], w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3627 OF 15495 ***
    // Wavefunction(s) for diagram number 3627
    // (none)
    // Amplitude(s) for diagram number 3627
    VVV1_0( w_fp[93], w_fp[100], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    VVV1_0( w_fp[380], w_fp[100], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    VVV1_0( w_fp[379], w_fp[100], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 3628 OF 15495 ***
    // Wavefunction(s) for diagram number 3628
    // (none)
    // Amplitude(s) for diagram number 3628
    FFV1_0( w_fp[461], w_fp[2], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[461], w_fp[2], w_fp[380], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[461], w_fp[2], w_fp[379], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3629 OF 15495 ***
    // Wavefunction(s) for diagram number 3629
    // (none)
    // Amplitude(s) for diagram number 3629
    VVV1_0( w_fp[479], w_fp[389], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[392], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[393], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];

    // *** DIAGRAM 3630 OF 15495 ***
    // Wavefunction(s) for diagram number 3630
    // (none)
    // Amplitude(s) for diagram number 3630
    VVV1_0( w_fp[479], w_fp[394], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[395], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[396], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 3631 OF 15495 ***
    // Wavefunction(s) for diagram number 3631
    // (none)
    // Amplitude(s) for diagram number 3631
    VVVV1_0( w_fp[97], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    VVVV3_0( w_fp[97], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    VVVV4_0( w_fp[97], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    VVVV1_0( w_fp[391], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    VVVV3_0( w_fp[391], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    VVVV4_0( w_fp[391], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    VVVV1_0( w_fp[390], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    VVVV3_0( w_fp[390], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    VVVV4_0( w_fp[390], w_fp[5], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];

    // *** DIAGRAM 3632 OF 15495 ***
    // Wavefunction(s) for diagram number 3632
    FFV1_1( w_fp[2], w_fp[97], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[461] );
    FFV1_1( w_fp[2], w_fp[391], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[469] );
    FFV1_1( w_fp[2], w_fp[390], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[487] );
    // Amplitude(s) for diagram number 3632
    FFV1_0( w_fp[442], w_fp[461], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[469], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[487], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];

    // *** DIAGRAM 3633 OF 15495 ***
    // Wavefunction(s) for diagram number 3633
    // (none)
    // Amplitude(s) for diagram number 3633
    FFV1_0( w_fp[442], w_fp[2], w_fp[394], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[395], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[396], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3634 OF 15495 ***
    // Wavefunction(s) for diagram number 3634
    // (none)
    // Amplitude(s) for diagram number 3634
    FFV1_0( w_fp[447], w_fp[461], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[469], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[487], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];

    // *** DIAGRAM 3635 OF 15495 ***
    // Wavefunction(s) for diagram number 3635
    // (none)
    // Amplitude(s) for diagram number 3635
    FFV1_0( w_fp[447], w_fp[2], w_fp[389], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[392], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[393], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3636 OF 15495 ***
    // Wavefunction(s) for diagram number 3636
    FFV1_2( w_fp[45], w_fp[97], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    FFV1_2( w_fp[45], w_fp[391], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[482] );
    FFV1_2( w_fp[45], w_fp[390], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[478] );
    // Amplitude(s) for diagram number 3636
    FFV1_0( w_fp[495], w_fp[169], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    FFV1_0( w_fp[482], w_fp[169], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    FFV1_0( w_fp[478], w_fp[169], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];

    // *** DIAGRAM 3637 OF 15495 ***
    // Wavefunction(s) for diagram number 3637
    // (none)
    // Amplitude(s) for diagram number 3637
    VVV1_0( w_fp[97], w_fp[6], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[391], w_fp[6], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[390], w_fp[6], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3638 OF 15495 ***
    // Wavefunction(s) for diagram number 3638
    // (none)
    // Amplitude(s) for diagram number 3638
    FFV1_0( w_fp[447], w_fp[169], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[169], w_fp[391], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[169], w_fp[390], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];

    // *** DIAGRAM 3639 OF 15495 ***
    // Wavefunction(s) for diagram number 3639
    // (none)
    // Amplitude(s) for diagram number 3639
    FFV1_0( w_fp[495], w_fp[197], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    FFV1_0( w_fp[482], w_fp[197], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    FFV1_0( w_fp[478], w_fp[197], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];

    // *** DIAGRAM 3640 OF 15495 ***
    // Wavefunction(s) for diagram number 3640
    // (none)
    // Amplitude(s) for diagram number 3640
    VVV1_0( w_fp[97], w_fp[5], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[391], w_fp[5], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[390], w_fp[5], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3641 OF 15495 ***
    // Wavefunction(s) for diagram number 3641
    // (none)
    // Amplitude(s) for diagram number 3641
    FFV1_0( w_fp[442], w_fp[197], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[197], w_fp[391], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[197], w_fp[390], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 3642 OF 15495 ***
    // Wavefunction(s) for diagram number 3642
    // (none)
    // Amplitude(s) for diagram number 3642
    FFV1_0( w_fp[495], w_fp[2], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[482], w_fp[2], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[478], w_fp[2], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3643 OF 15495 ***
    // Wavefunction(s) for diagram number 3643
    // (none)
    // Amplitude(s) for diagram number 3643
    VVV1_0( w_fp[97], w_fp[113], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    VVV1_0( w_fp[391], w_fp[113], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    VVV1_0( w_fp[390], w_fp[113], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];

    // *** DIAGRAM 3644 OF 15495 ***
    // Wavefunction(s) for diagram number 3644
    // (none)
    // Amplitude(s) for diagram number 3644
    FFV1_0( w_fp[464], w_fp[2], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[464], w_fp[2], w_fp[391], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[464], w_fp[2], w_fp[390], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3645 OF 15495 ***
    // Wavefunction(s) for diagram number 3645
    // (none)
    // Amplitude(s) for diagram number 3645
    VVV1_0( w_fp[479], w_fp[400], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[403], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[404], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3646 OF 15495 ***
    // Wavefunction(s) for diagram number 3646
    // (none)
    // Amplitude(s) for diagram number 3646
    VVV1_0( w_fp[479], w_fp[405], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[406], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[407], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3647 OF 15495 ***
    // Wavefunction(s) for diagram number 3647
    // (none)
    // Amplitude(s) for diagram number 3647
    VVVV1_0( w_fp[127], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVVV3_0( w_fp[127], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVVV4_0( w_fp[127], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    VVVV1_0( w_fp[402], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV3_0( w_fp[402], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV4_0( w_fp[402], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    VVVV1_0( w_fp[401], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV3_0( w_fp[401], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0( w_fp[401], w_fp[4], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 3648 OF 15495 ***
    // Wavefunction(s) for diagram number 3648
    FFV1_1( w_fp[2], w_fp[127], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[464] );
    FFV1_1( w_fp[2], w_fp[402], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[478] );
    FFV1_1( w_fp[2], w_fp[401], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[482] );
    // Amplitude(s) for diagram number 3648
    FFV1_0( w_fp[435], w_fp[464], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[478], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[482], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];

    // *** DIAGRAM 3649 OF 15495 ***
    // Wavefunction(s) for diagram number 3649
    // (none)
    // Amplitude(s) for diagram number 3649
    FFV1_0( w_fp[435], w_fp[2], w_fp[405], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[406], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[407], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3650 OF 15495 ***
    // Wavefunction(s) for diagram number 3650
    // (none)
    // Amplitude(s) for diagram number 3650
    FFV1_0( w_fp[451], w_fp[464], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[478], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[482], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[561] -= amp_sv[0];

    // *** DIAGRAM 3651 OF 15495 ***
    // Wavefunction(s) for diagram number 3651
    // (none)
    // Amplitude(s) for diagram number 3651
    FFV1_0( w_fp[451], w_fp[2], w_fp[400], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[403], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[404], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3652 OF 15495 ***
    // Wavefunction(s) for diagram number 3652
    FFV1_2( w_fp[45], w_fp[127], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    FFV1_2( w_fp[45], w_fp[402], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[470] );
    FFV1_2( w_fp[45], w_fp[401], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[511] );
    // Amplitude(s) for diagram number 3652
    FFV1_0( w_fp[495], w_fp[156], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    FFV1_0( w_fp[470], w_fp[156], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    FFV1_0( w_fp[511], w_fp[156], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];

    // *** DIAGRAM 3653 OF 15495 ***
    // Wavefunction(s) for diagram number 3653
    // (none)
    // Amplitude(s) for diagram number 3653
    VVV1_0( w_fp[127], w_fp[7], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[402], w_fp[7], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[279] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[401], w_fp[7], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3654 OF 15495 ***
    // Wavefunction(s) for diagram number 3654
    // (none)
    // Amplitude(s) for diagram number 3654
    FFV1_0( w_fp[451], w_fp[156], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] += amp_sv[0];
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[156], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[279] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[156], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];

    // *** DIAGRAM 3655 OF 15495 ***
    // Wavefunction(s) for diagram number 3655
    // (none)
    // Amplitude(s) for diagram number 3655
    FFV1_0( w_fp[495], w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    FFV1_0( w_fp[470], w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    FFV1_0( w_fp[511], w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 3656 OF 15495 ***
    // Wavefunction(s) for diagram number 3656
    // (none)
    // Amplitude(s) for diagram number 3656
    VVV1_0( w_fp[127], w_fp[4], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[402], w_fp[4], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[401], w_fp[4], w_fp[492], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3657 OF 15495 ***
    // Wavefunction(s) for diagram number 3657
    // (none)
    // Amplitude(s) for diagram number 3657
    FFV1_0( w_fp[435], w_fp[215], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[215], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[215], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3658 OF 15495 ***
    // Wavefunction(s) for diagram number 3658
    // (none)
    // Amplitude(s) for diagram number 3658
    FFV1_0( w_fp[495], w_fp[2], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[470], w_fp[2], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[511], w_fp[2], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3659 OF 15495 ***
    // Wavefunction(s) for diagram number 3659
    // (none)
    // Amplitude(s) for diagram number 3659
    VVV1_0( w_fp[127], w_fp[102], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[561] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    VVV1_0( w_fp[402], w_fp[102], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[347] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    VVV1_0( w_fp[401], w_fp[102], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[345] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 3660 OF 15495 ***
    // Wavefunction(s) for diagram number 3660
    // (none)
    // Amplitude(s) for diagram number 3660
    FFV1_0( w_fp[462], w_fp[2], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[462], w_fp[2], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[462], w_fp[2], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3661 OF 15495 ***
    // Wavefunction(s) for diagram number 3661
    // (none)
    // Amplitude(s) for diagram number 3661
    VVV1_0( w_fp[479], w_fp[411], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[284], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[414], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];

    // *** DIAGRAM 3662 OF 15495 ***
    // Wavefunction(s) for diagram number 3662
    // (none)
    // Amplitude(s) for diagram number 3662
    VVV1_0( w_fp[479], w_fp[415], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[416], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[417], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 3663 OF 15495 ***
    // Wavefunction(s) for diagram number 3663
    // (none)
    // Amplitude(s) for diagram number 3663
    VVVV1_0( w_fp[119], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    VVVV3_0( w_fp[119], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    VVVV4_0( w_fp[119], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    VVVV1_0( w_fp[413], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    VVVV3_0( w_fp[413], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    VVVV4_0( w_fp[413], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    VVVV1_0( w_fp[412], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    VVVV3_0( w_fp[412], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVVV4_0( w_fp[412], w_fp[4], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 3664 OF 15495 ***
    // Wavefunction(s) for diagram number 3664
    FFV1_1( w_fp[2], w_fp[119], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[462] );
    FFV1_1( w_fp[2], w_fp[413], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[511] );
    FFV1_1( w_fp[2], w_fp[412], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[470] );
    // Amplitude(s) for diagram number 3664
    FFV1_0( w_fp[435], w_fp[462], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[511], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[470], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];

    // *** DIAGRAM 3665 OF 15495 ***
    // Wavefunction(s) for diagram number 3665
    // (none)
    // Amplitude(s) for diagram number 3665
    FFV1_0( w_fp[435], w_fp[2], w_fp[415], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[416], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[417], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3666 OF 15495 ***
    // Wavefunction(s) for diagram number 3666
    // (none)
    // Amplitude(s) for diagram number 3666
    FFV1_0( w_fp[447], w_fp[462], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[511], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[470], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];

    // *** DIAGRAM 3667 OF 15495 ***
    // Wavefunction(s) for diagram number 3667
    // (none)
    // Amplitude(s) for diagram number 3667
    FFV1_0( w_fp[447], w_fp[2], w_fp[411], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[284], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[414], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3668 OF 15495 ***
    // Wavefunction(s) for diagram number 3668
    FFV1_2( w_fp[45], w_fp[119], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    FFV1_2( w_fp[45], w_fp[413], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[492] );
    FFV1_2( w_fp[45], w_fp[412], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[485] );
    // Amplitude(s) for diagram number 3668
    FFV1_0( w_fp[495], w_fp[156], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    FFV1_0( w_fp[492], w_fp[156], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    FFV1_0( w_fp[485], w_fp[156], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];

    // *** DIAGRAM 3669 OF 15495 ***
    // Wavefunction(s) for diagram number 3669
    // (none)
    // Amplitude(s) for diagram number 3669
    VVV1_0( w_fp[119], w_fp[6], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[413], w_fp[6], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[285] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[412], w_fp[6], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[345] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3670 OF 15495 ***
    // Wavefunction(s) for diagram number 3670
    // (none)
    // Amplitude(s) for diagram number 3670
    FFV1_0( w_fp[447], w_fp[156], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] += amp_sv[0];
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[351] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[156], w_fp[413], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[285] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[156], w_fp[412], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];

    // *** DIAGRAM 3671 OF 15495 ***
    // Wavefunction(s) for diagram number 3671
    // (none)
    // Amplitude(s) for diagram number 3671
    FFV1_0( w_fp[495], w_fp[197], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    FFV1_0( w_fp[492], w_fp[197], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    FFV1_0( w_fp[485], w_fp[197], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];

    // *** DIAGRAM 3672 OF 15495 ***
    // Wavefunction(s) for diagram number 3672
    // (none)
    // Amplitude(s) for diagram number 3672
    VVV1_0( w_fp[119], w_fp[4], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[413], w_fp[4], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[412], w_fp[4], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3673 OF 15495 ***
    // Wavefunction(s) for diagram number 3673
    // (none)
    // Amplitude(s) for diagram number 3673
    FFV1_0( w_fp[435], w_fp[197], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[197], w_fp[413], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[197], w_fp[412], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];

    // *** DIAGRAM 3674 OF 15495 ***
    // Wavefunction(s) for diagram number 3674
    // (none)
    // Amplitude(s) for diagram number 3674
    FFV1_0( w_fp[495], w_fp[2], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[492], w_fp[2], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[485], w_fp[2], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[321] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3675 OF 15495 ***
    // Wavefunction(s) for diagram number 3675
    // (none)
    // Amplitude(s) for diagram number 3675
    VVV1_0( w_fp[119], w_fp[86], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    VVV1_0( w_fp[413], w_fp[86], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[323] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    VVV1_0( w_fp[412], w_fp[86], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[321] -= amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];

    // *** DIAGRAM 3676 OF 15495 ***
    // Wavefunction(s) for diagram number 3676
    // (none)
    // Amplitude(s) for diagram number 3676
    FFV1_0( w_fp[459], w_fp[2], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[459], w_fp[2], w_fp[413], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[459], w_fp[2], w_fp[412], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3677 OF 15495 ***
    // Wavefunction(s) for diagram number 3677
    // (none)
    // Amplitude(s) for diagram number 3677
    VVV1_0( w_fp[479], w_fp[421], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[282], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[424], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 3678 OF 15495 ***
    // Wavefunction(s) for diagram number 3678
    // (none)
    // Amplitude(s) for diagram number 3678
    VVV1_0( w_fp[479], w_fp[425], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[426], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[239] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[427], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 3679 OF 15495 ***
    // Wavefunction(s) for diagram number 3679
    // (none)
    // Amplitude(s) for diagram number 3679
    VVVV1_0( w_fp[111], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    VVVV3_0( w_fp[111], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    VVVV4_0( w_fp[111], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    VVVV1_0( w_fp[423], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    VVVV3_0( w_fp[423], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[239] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    VVVV4_0( w_fp[423], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    VVVV1_0( w_fp[422], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    VVVV3_0( w_fp[422], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVVV4_0( w_fp[422], w_fp[4], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 3680 OF 15495 ***
    // Wavefunction(s) for diagram number 3680
    FFV1_1( w_fp[2], w_fp[111], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[459] );
    FFV1_1( w_fp[2], w_fp[423], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[485] );
    FFV1_1( w_fp[2], w_fp[422], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[492] );
    // Amplitude(s) for diagram number 3680
    FFV1_0( w_fp[435], w_fp[459], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[485], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[492], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

    // *** DIAGRAM 3681 OF 15495 ***
    // Wavefunction(s) for diagram number 3681
    // (none)
    // Amplitude(s) for diagram number 3681
    FFV1_0( w_fp[435], w_fp[2], w_fp[425], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[426], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[427], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3682 OF 15495 ***
    // Wavefunction(s) for diagram number 3682
    // (none)
    // Amplitude(s) for diagram number 3682
    FFV1_0( w_fp[442], w_fp[459], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[485], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[492], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];

    // *** DIAGRAM 3683 OF 15495 ***
    // Wavefunction(s) for diagram number 3683
    // (none)
    // Amplitude(s) for diagram number 3683
    FFV1_0( w_fp[442], w_fp[2], w_fp[421], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[282], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[424], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3684 OF 15495 ***
    // Wavefunction(s) for diagram number 3684
    FFV1_2( w_fp[45], w_fp[111], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    FFV1_2( w_fp[45], w_fp[423], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[490] );
    FFV1_2( w_fp[45], w_fp[422], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[514] );
    // Amplitude(s) for diagram number 3684
    FFV1_0( w_fp[495], w_fp[156], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    FFV1_0( w_fp[490], w_fp[156], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    FFV1_0( w_fp[514], w_fp[156], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];

    // *** DIAGRAM 3685 OF 15495 ***
    // Wavefunction(s) for diagram number 3685
    // (none)
    // Amplitude(s) for diagram number 3685
    VVV1_0( w_fp[111], w_fp[5], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[423], w_fp[5], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[287] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[422], w_fp[5], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[323] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[347] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3686 OF 15495 ***
    // Wavefunction(s) for diagram number 3686
    // (none)
    // Amplitude(s) for diagram number 3686
    FFV1_0( w_fp[442], w_fp[156], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[156], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[156], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[323] += amp_sv[0];
    jamp_sv[347] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];

    // *** DIAGRAM 3687 OF 15495 ***
    // Wavefunction(s) for diagram number 3687
    // (none)
    // Amplitude(s) for diagram number 3687
    FFV1_0( w_fp[495], w_fp[169], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    FFV1_0( w_fp[490], w_fp[169], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    FFV1_0( w_fp[514], w_fp[169], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];

    // *** DIAGRAM 3688 OF 15495 ***
    // Wavefunction(s) for diagram number 3688
    // (none)
    // Amplitude(s) for diagram number 3688
    VVV1_0( w_fp[111], w_fp[4], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[423], w_fp[4], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[422], w_fp[4], w_fp[488], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3689 OF 15495 ***
    // Wavefunction(s) for diagram number 3689
    // (none)
    // Amplitude(s) for diagram number 3689
    FFV1_0( w_fp[435], w_fp[169], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[169], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[169], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];

    // *** DIAGRAM 3690 OF 15495 ***
    // Wavefunction(s) for diagram number 3690
    // (none)
    // Amplitude(s) for diagram number 3690
    FFV1_0( w_fp[495], w_fp[2], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[490], w_fp[2], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[299] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[514], w_fp[2], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[297] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3691 OF 15495 ***
    // Wavefunction(s) for diagram number 3691
    // (none)
    // Amplitude(s) for diagram number 3691
    VVV1_0( w_fp[111], w_fp[66], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVV1_0( w_fp[423], w_fp[66], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    VVV1_0( w_fp[422], w_fp[66], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[297] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

    // *** DIAGRAM 3692 OF 15495 ***
    // Wavefunction(s) for diagram number 3692
    // (none)
    // Amplitude(s) for diagram number 3692
    FFV1_0( w_fp[456], w_fp[2], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[456], w_fp[2], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[456], w_fp[2], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3693 OF 15495 ***
    // Wavefunction(s) for diagram number 3693
    // (none)
    // Amplitude(s) for diagram number 3693
    FFV1_0( w_fp[505], w_fp[195], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[67], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[15], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];

    // *** DIAGRAM 3694 OF 15495 ***
    // Wavefunction(s) for diagram number 3694
    // (none)
    // Amplitude(s) for diagram number 3694
    FFV1_0( w_fp[505], w_fp[2], w_fp[165], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[240], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3695 OF 15495 ***
    // Wavefunction(s) for diagram number 3695
    // (none)
    // Amplitude(s) for diagram number 3695
    VVV1_0( w_fp[479], w_fp[91], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[433], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[432], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3696 OF 15495 ***
    // Wavefunction(s) for diagram number 3696
    // (none)
    // Amplitude(s) for diagram number 3696
    VVV1_0( w_fp[479], w_fp[1], w_fp[165], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[240], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[30], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3697 OF 15495 ***
    // Wavefunction(s) for diagram number 3697
    // (none)
    // Amplitude(s) for diagram number 3697
    VVVV1_0( w_fp[1], w_fp[133], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[133], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[133], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[134], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[134], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[134], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[135], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[135], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[135], w_fp[7], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3698 OF 15495 ***
    // Wavefunction(s) for diagram number 3698
    // (none)
    // Amplitude(s) for diagram number 3698
    FFV1_0( w_fp[451], w_fp[2], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[433], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[451], w_fp[2], w_fp[432], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3699 OF 15495 ***
    // Wavefunction(s) for diagram number 3699
    // (none)
    // Amplitude(s) for diagram number 3699
    FFV1_0( w_fp[451], w_fp[195], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[303] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[67], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[15], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];

    // *** DIAGRAM 3700 OF 15495 ***
    // Wavefunction(s) for diagram number 3700
    // (none)
    // Amplitude(s) for diagram number 3700
    FFV1_0( w_fp[505], w_fp[51], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[166], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[241], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 3701 OF 15495 ***
    // Wavefunction(s) for diagram number 3701
    // (none)
    // Amplitude(s) for diagram number 3701
    FFV1_0( w_fp[505], w_fp[2], w_fp[162], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3702 OF 15495 ***
    // Wavefunction(s) for diagram number 3702
    // (none)
    // Amplitude(s) for diagram number 3702
    VVV1_0( w_fp[479], w_fp[108], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[431], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[280], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 3703 OF 15495 ***
    // Wavefunction(s) for diagram number 3703
    // (none)
    // Amplitude(s) for diagram number 3703
    VVV1_0( w_fp[479], w_fp[1], w_fp[162], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[53], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[52], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 3704 OF 15495 ***
    // Wavefunction(s) for diagram number 3704
    // (none)
    // Amplitude(s) for diagram number 3704
    VVVV1_0( w_fp[1], w_fp[139], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[139], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[139], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[140], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[140], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[140], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[141], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[141], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[141], w_fp[6], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];

    // *** DIAGRAM 3705 OF 15495 ***
    // Wavefunction(s) for diagram number 3705
    // (none)
    // Amplitude(s) for diagram number 3705
    FFV1_0( w_fp[447], w_fp[2], w_fp[108], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[431], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[447], w_fp[2], w_fp[280], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3706 OF 15495 ***
    // Wavefunction(s) for diagram number 3706
    // (none)
    // Amplitude(s) for diagram number 3706
    FFV1_0( w_fp[447], w_fp[51], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[309] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[166], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[241], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];

    // *** DIAGRAM 3707 OF 15495 ***
    // Wavefunction(s) for diagram number 3707
    // (none)
    // Amplitude(s) for diagram number 3707
    FFV1_0( w_fp[505], w_fp[242], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[47], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[13], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 3708 OF 15495 ***
    // Wavefunction(s) for diagram number 3708
    // (none)
    // Amplitude(s) for diagram number 3708
    FFV1_0( w_fp[505], w_fp[2], w_fp[159], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[167], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3709 OF 15495 ***
    // Wavefunction(s) for diagram number 3709
    // (none)
    // Amplitude(s) for diagram number 3709
    VVV1_0( w_fp[479], w_fp[177], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[305], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[10], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3710 OF 15495 ***
    // Wavefunction(s) for diagram number 3710
    // (none)
    // Amplitude(s) for diagram number 3710
    VVV1_0( w_fp[479], w_fp[1], w_fp[159], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[54], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[167], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 3711 OF 15495 ***
    // Wavefunction(s) for diagram number 3711
    // (none)
    // Amplitude(s) for diagram number 3711
    VVVV1_0( w_fp[1], w_fp[145], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[145], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[145], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[146], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[146], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[146], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[147], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[147], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[147], w_fp[5], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 3712 OF 15495 ***
    // Wavefunction(s) for diagram number 3712
    // (none)
    // Amplitude(s) for diagram number 3712
    FFV1_0( w_fp[442], w_fp[2], w_fp[177], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[305], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[442], w_fp[2], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3713 OF 15495 ***
    // Wavefunction(s) for diagram number 3713
    // (none)
    // Amplitude(s) for diagram number 3713
    FFV1_0( w_fp[442], w_fp[242], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] += amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[47], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[13], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 3714 OF 15495 ***
    // Wavefunction(s) for diagram number 3714
    // (none)
    // Amplitude(s) for diagram number 3714
    FFV1_0( w_fp[505], w_fp[229], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[28], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    FFV1_0( w_fp[505], w_fp[60], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3715 OF 15495 ***
    // Wavefunction(s) for diagram number 3715
    // (none)
    // Amplitude(s) for diagram number 3715
    FFV1_0( w_fp[505], w_fp[2], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[505], w_fp[2], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3716 OF 15495 ***
    // Wavefunction(s) for diagram number 3716
    // (none)
    // Amplitude(s) for diagram number 3716
    VVV1_0( w_fp[479], w_fp[345], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[95], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    VVV1_0( w_fp[479], w_fp[434], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3717 OF 15495 ***
    // Wavefunction(s) for diagram number 3717
    // (none)
    // Amplitude(s) for diagram number 3717
    VVV1_0( w_fp[479], w_fp[1], w_fp[58], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[57], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    VVV1_0( w_fp[479], w_fp[1], w_fp[38], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 3718 OF 15495 ***
    // Wavefunction(s) for diagram number 3718
    // (none)
    // Amplitude(s) for diagram number 3718
    VVVV1_0( w_fp[1], w_fp[4], w_fp[151], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[4], w_fp[151], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[273] += amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[287] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[4], w_fp[151], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[4], w_fp[152], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[4], w_fp[152], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[275] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[281] -= amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[4], w_fp[152], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[4], w_fp[153], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[4], w_fp[153], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[279] += amp_sv[0];
    jamp_sv[285] += amp_sv[0];
    jamp_sv[287] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[4], w_fp[153], w_fp[479], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[159] += amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3719 OF 15495 ***
    // Wavefunction(s) for diagram number 3719
    // (none)
    // Amplitude(s) for diagram number 3719
    FFV1_0( w_fp[435], w_fp[2], w_fp[345], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[95], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[435], w_fp[2], w_fp[434], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3720 OF 15495 ***
    // Wavefunction(s) for diagram number 3720
    // (none)
    // Amplitude(s) for diagram number 3720
    FFV1_0( w_fp[435], w_fp[229], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[28], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    FFV1_0( w_fp[435], w_fp[60], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 3721 OF 15495 ***
    // Wavefunction(s) for diagram number 3721
    VVV1P0_1( w_fp[0], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[435] );
    FFV1_2( w_fp[3], w_fp[435], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[479] );
    FFV1_2( w_fp[479], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[505] );
    // Amplitude(s) for diagram number 3721
    FFV1_0( w_fp[505], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3722 OF 15495 ***
    // Wavefunction(s) for diagram number 3722
    FFV1_2( w_fp[479], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[442] );
    // Amplitude(s) for diagram number 3722
    FFV1_0( w_fp[442], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3723 OF 15495 ***
    // Wavefunction(s) for diagram number 3723
    FFV1_2( w_fp[479], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[447] );
    // Amplitude(s) for diagram number 3723
    FFV1_0( w_fp[447], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3724 OF 15495 ***
    // Wavefunction(s) for diagram number 3724
    // (none)
    // Amplitude(s) for diagram number 3724
    FFV1_0( w_fp[442], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3725 OF 15495 ***
    // Wavefunction(s) for diagram number 3725
    // (none)
    // Amplitude(s) for diagram number 3725
    FFV1_0( w_fp[447], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3726 OF 15495 ***
    // Wavefunction(s) for diagram number 3726
    // (none)
    // Amplitude(s) for diagram number 3726
    FFV1_0( w_fp[505], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3727 OF 15495 ***
    // Wavefunction(s) for diagram number 3727
    VVV1P0_1( w_fp[435], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[451] );
    FFV1P0_3( w_fp[3], w_fp[259], COUPs[1], 1.0, depCoup, 0., 0., w_fp[456] );
    // Amplitude(s) for diagram number 3727
    VVVV1_0( w_fp[451], w_fp[456], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVVV3_0( w_fp[451], w_fp[456], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVVV4_0( w_fp[451], w_fp[456], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 3728 OF 15495 ***
    // Wavefunction(s) for diagram number 3728
    VVV1P0_1( w_fp[451], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[514] );
    // Amplitude(s) for diagram number 3728
    VVV1_0( w_fp[456], w_fp[7], w_fp[514], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 3729 OF 15495 ***
    // Wavefunction(s) for diagram number 3729
    VVV1P0_1( w_fp[451], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[490] );
    // Amplitude(s) for diagram number 3729
    VVV1_0( w_fp[456], w_fp[6], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 3730 OF 15495 ***
    // Wavefunction(s) for diagram number 3730
    FFV1_2( w_fp[3], w_fp[451], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    // Amplitude(s) for diagram number 3730
    FFV1_0( w_fp[495], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];

    // *** DIAGRAM 3731 OF 15495 ***
    // Wavefunction(s) for diagram number 3731
    // (none)
    // Amplitude(s) for diagram number 3731
    FFV1_0( w_fp[3], w_fp[439], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3732 OF 15495 ***
    // Wavefunction(s) for diagram number 3732
    // (none)
    // Amplitude(s) for diagram number 3732
    FFV1_0( w_fp[495], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 3733 OF 15495 ***
    // Wavefunction(s) for diagram number 3733
    // (none)
    // Amplitude(s) for diagram number 3733
    FFV1_0( w_fp[3], w_fp[441], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3734 OF 15495 ***
    // Wavefunction(s) for diagram number 3734
    VVV1P0_1( w_fp[435], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[488] );
    // Amplitude(s) for diagram number 3734
    VVVV1_0( w_fp[488], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 3735 OF 15495 ***
    // Wavefunction(s) for diagram number 3735
    VVV1P0_1( w_fp[488], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[497] );
    // Amplitude(s) for diagram number 3735
    VVV1_0( w_fp[456], w_fp[7], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 3736 OF 15495 ***
    // Wavefunction(s) for diagram number 3736
    VVV1P0_1( w_fp[488], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[45] );
    // Amplitude(s) for diagram number 3736
    VVV1_0( w_fp[456], w_fp[5], w_fp[45], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 3737 OF 15495 ***
    // Wavefunction(s) for diagram number 3737
    FFV1_2( w_fp[3], w_fp[488], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[450] );
    // Amplitude(s) for diagram number 3737
    FFV1_0( w_fp[450], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];

    // *** DIAGRAM 3738 OF 15495 ***
    // Wavefunction(s) for diagram number 3738
    // (none)
    // Amplitude(s) for diagram number 3738
    FFV1_0( w_fp[3], w_fp[436], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3739 OF 15495 ***
    // Wavefunction(s) for diagram number 3739
    // (none)
    // Amplitude(s) for diagram number 3739
    FFV1_0( w_fp[450], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];

    // *** DIAGRAM 3740 OF 15495 ***
    // Wavefunction(s) for diagram number 3740
    // (none)
    // Amplitude(s) for diagram number 3740
    FFV1_0( w_fp[3], w_fp[441], w_fp[497], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3741 OF 15495 ***
    // Wavefunction(s) for diagram number 3741
    VVV1P0_1( w_fp[435], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[437] );
    // Amplitude(s) for diagram number 3741
    VVVV1_0( w_fp[437], w_fp[456], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    VVVV3_0( w_fp[437], w_fp[456], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    VVVV4_0( w_fp[437], w_fp[456], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];

    // *** DIAGRAM 3742 OF 15495 ***
    // Wavefunction(s) for diagram number 3742
    VVV1P0_1( w_fp[437], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[480] );
    // Amplitude(s) for diagram number 3742
    VVV1_0( w_fp[456], w_fp[6], w_fp[480], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];

    // *** DIAGRAM 3743 OF 15495 ***
    // Wavefunction(s) for diagram number 3743
    VVV1P0_1( w_fp[437], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[509] );
    // Amplitude(s) for diagram number 3743
    VVV1_0( w_fp[456], w_fp[5], w_fp[509], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];

    // *** DIAGRAM 3744 OF 15495 ***
    // Wavefunction(s) for diagram number 3744
    FFV1_2( w_fp[3], w_fp[437], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[486] );
    // Amplitude(s) for diagram number 3744
    FFV1_0( w_fp[486], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];

    // *** DIAGRAM 3745 OF 15495 ***
    // Wavefunction(s) for diagram number 3745
    // (none)
    // Amplitude(s) for diagram number 3745
    FFV1_0( w_fp[3], w_fp[436], w_fp[509], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3746 OF 15495 ***
    // Wavefunction(s) for diagram number 3746
    // (none)
    // Amplitude(s) for diagram number 3746
    FFV1_0( w_fp[486], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];

    // *** DIAGRAM 3747 OF 15495 ***
    // Wavefunction(s) for diagram number 3747
    // (none)
    // Amplitude(s) for diagram number 3747
    FFV1_0( w_fp[3], w_fp[439], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3748 OF 15495 ***
    // Wavefunction(s) for diagram number 3748
    VVVV1P0_1( w_fp[435], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[515] );
    VVVV3P0_1( w_fp[435], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[454] );
    VVVV4P0_1( w_fp[435], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[438] );
    // Amplitude(s) for diagram number 3748
    VVV1_0( w_fp[515], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVV1_0( w_fp[454], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    VVV1_0( w_fp[438], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 3749 OF 15495 ***
    // Wavefunction(s) for diagram number 3749
    // (none)
    // Amplitude(s) for diagram number 3749
    FFV1_0( w_fp[3], w_fp[441], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3750 OF 15495 ***
    // Wavefunction(s) for diagram number 3750
    VVVV1P0_1( w_fp[435], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[481] );
    VVVV3P0_1( w_fp[435], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[510] );
    VVVV4P0_1( w_fp[435], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[448] );
    // Amplitude(s) for diagram number 3750
    VVV1_0( w_fp[481], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    VVV1_0( w_fp[510], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    VVV1_0( w_fp[448], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 3751 OF 15495 ***
    // Wavefunction(s) for diagram number 3751
    // (none)
    // Amplitude(s) for diagram number 3751
    FFV1_0( w_fp[3], w_fp[439], w_fp[481], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[439], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[439], w_fp[448], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3752 OF 15495 ***
    // Wavefunction(s) for diagram number 3752
    VVVV1P0_1( w_fp[435], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[444] );
    VVVV3P0_1( w_fp[435], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[452] );
    VVVV4P0_1( w_fp[435], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[445] );
    // Amplitude(s) for diagram number 3752
    VVV1_0( w_fp[444], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    VVV1_0( w_fp[452], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    VVV1_0( w_fp[445], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 3753 OF 15495 ***
    // Wavefunction(s) for diagram number 3753
    // (none)
    // Amplitude(s) for diagram number 3753
    FFV1_0( w_fp[3], w_fp[436], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[436], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[436], w_fp[445], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3754 OF 15495 ***
    // Wavefunction(s) for diagram number 3754
    FFV1_1( w_fp[259], w_fp[435], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[453] );
    // Amplitude(s) for diagram number 3754
    FFV1_0( w_fp[170], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3755 OF 15495 ***
    // Wavefunction(s) for diagram number 3755
    // (none)
    // Amplitude(s) for diagram number 3755
    FFV1_0( w_fp[171], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3756 OF 15495 ***
    // Wavefunction(s) for diagram number 3756
    FFV1_2( w_fp[168], w_fp[435], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[449] );
    // Amplitude(s) for diagram number 3756
    FFV1_0( w_fp[449], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3757 OF 15495 ***
    // Wavefunction(s) for diagram number 3757
    // (none)
    // Amplitude(s) for diagram number 3757
    FFV1_0( w_fp[449], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3758 OF 15495 ***
    // Wavefunction(s) for diagram number 3758
    FFV1P0_3( w_fp[168], w_fp[259], COUPs[1], 1.0, depCoup, 0., 0., w_fp[468] );
    // Amplitude(s) for diagram number 3758
    VVV1_0( w_fp[488], w_fp[468], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3759 OF 15495 ***
    // Wavefunction(s) for diagram number 3759
    // (none)
    // Amplitude(s) for diagram number 3759
    FFV1_0( w_fp[168], w_fp[441], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 3760 OF 15495 ***
    // Wavefunction(s) for diagram number 3760
    // (none)
    // Amplitude(s) for diagram number 3760
    FFV1_0( w_fp[171], w_fp[259], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];

    // *** DIAGRAM 3761 OF 15495 ***
    // Wavefunction(s) for diagram number 3761
    // (none)
    // Amplitude(s) for diagram number 3761
    VVV1_0( w_fp[437], w_fp[468], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3762 OF 15495 ***
    // Wavefunction(s) for diagram number 3762
    // (none)
    // Amplitude(s) for diagram number 3762
    FFV1_0( w_fp[168], w_fp[439], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];

    // *** DIAGRAM 3763 OF 15495 ***
    // Wavefunction(s) for diagram number 3763
    // (none)
    // Amplitude(s) for diagram number 3763
    FFV1_0( w_fp[170], w_fp[259], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];

    // *** DIAGRAM 3764 OF 15495 ***
    // Wavefunction(s) for diagram number 3764
    // (none)
    // Amplitude(s) for diagram number 3764
    FFV1_0( w_fp[171], w_fp[439], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3765 OF 15495 ***
    // Wavefunction(s) for diagram number 3765
    // (none)
    // Amplitude(s) for diagram number 3765
    FFV1_0( w_fp[170], w_fp[441], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3766 OF 15495 ***
    // Wavefunction(s) for diagram number 3766
    // (none)
    // Amplitude(s) for diagram number 3766
    FFV1_0( w_fp[168], w_fp[259], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[445], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3767 OF 15495 ***
    // Wavefunction(s) for diagram number 3767
    // (none)
    // Amplitude(s) for diagram number 3767
    FFV1_0( w_fp[168], w_fp[453], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];

    // *** DIAGRAM 3768 OF 15495 ***
    // Wavefunction(s) for diagram number 3768
    // (none)
    // Amplitude(s) for diagram number 3768
    FFV1_0( w_fp[449], w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 3769 OF 15495 ***
    // Wavefunction(s) for diagram number 3769
    VVV1P0_1( w_fp[435], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[472] );
    // Amplitude(s) for diagram number 3769
    FFV1_0( w_fp[168], w_fp[259], w_fp[472], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3770 OF 15495 ***
    // Wavefunction(s) for diagram number 3770
    // (none)
    // Amplitude(s) for diagram number 3770
    FFV1_0( w_fp[175], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3771 OF 15495 ***
    // Wavefunction(s) for diagram number 3771
    // (none)
    // Amplitude(s) for diagram number 3771
    FFV1_0( w_fp[176], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3772 OF 15495 ***
    // Wavefunction(s) for diagram number 3772
    FFV1_2( w_fp[174], w_fp[435], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[471] );
    // Amplitude(s) for diagram number 3772
    FFV1_0( w_fp[471], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3773 OF 15495 ***
    // Wavefunction(s) for diagram number 3773
    // (none)
    // Amplitude(s) for diagram number 3773
    FFV1_0( w_fp[471], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3774 OF 15495 ***
    // Wavefunction(s) for diagram number 3774
    FFV1P0_3( w_fp[174], w_fp[259], COUPs[1], 1.0, depCoup, 0., 0., w_fp[475] );
    // Amplitude(s) for diagram number 3774
    VVV1_0( w_fp[451], w_fp[475], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3775 OF 15495 ***
    // Wavefunction(s) for diagram number 3775
    // (none)
    // Amplitude(s) for diagram number 3775
    FFV1_0( w_fp[174], w_fp[441], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 3776 OF 15495 ***
    // Wavefunction(s) for diagram number 3776
    // (none)
    // Amplitude(s) for diagram number 3776
    FFV1_0( w_fp[176], w_fp[259], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];

    // *** DIAGRAM 3777 OF 15495 ***
    // Wavefunction(s) for diagram number 3777
    // (none)
    // Amplitude(s) for diagram number 3777
    VVV1_0( w_fp[437], w_fp[475], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3778 OF 15495 ***
    // Wavefunction(s) for diagram number 3778
    // (none)
    // Amplitude(s) for diagram number 3778
    FFV1_0( w_fp[174], w_fp[436], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];

    // *** DIAGRAM 3779 OF 15495 ***
    // Wavefunction(s) for diagram number 3779
    // (none)
    // Amplitude(s) for diagram number 3779
    FFV1_0( w_fp[175], w_fp[259], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];

    // *** DIAGRAM 3780 OF 15495 ***
    // Wavefunction(s) for diagram number 3780
    // (none)
    // Amplitude(s) for diagram number 3780
    FFV1_0( w_fp[176], w_fp[436], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3781 OF 15495 ***
    // Wavefunction(s) for diagram number 3781
    // (none)
    // Amplitude(s) for diagram number 3781
    FFV1_0( w_fp[175], w_fp[441], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3782 OF 15495 ***
    // Wavefunction(s) for diagram number 3782
    // (none)
    // Amplitude(s) for diagram number 3782
    FFV1_0( w_fp[174], w_fp[259], w_fp[481], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[448], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3783 OF 15495 ***
    // Wavefunction(s) for diagram number 3783
    // (none)
    // Amplitude(s) for diagram number 3783
    FFV1_0( w_fp[174], w_fp[453], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];

    // *** DIAGRAM 3784 OF 15495 ***
    // Wavefunction(s) for diagram number 3784
    // (none)
    // Amplitude(s) for diagram number 3784
    FFV1_0( w_fp[471], w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 3785 OF 15495 ***
    // Wavefunction(s) for diagram number 3785
    VVV1P0_1( w_fp[435], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[474] );
    // Amplitude(s) for diagram number 3785
    FFV1_0( w_fp[174], w_fp[259], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3786 OF 15495 ***
    // Wavefunction(s) for diagram number 3786
    // (none)
    // Amplitude(s) for diagram number 3786
    FFV1_0( w_fp[180], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3787 OF 15495 ***
    // Wavefunction(s) for diagram number 3787
    // (none)
    // Amplitude(s) for diagram number 3787
    FFV1_0( w_fp[181], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3788 OF 15495 ***
    // Wavefunction(s) for diagram number 3788
    FFV1_2( w_fp[179], w_fp[435], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[516] );
    // Amplitude(s) for diagram number 3788
    FFV1_0( w_fp[516], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3789 OF 15495 ***
    // Wavefunction(s) for diagram number 3789
    // (none)
    // Amplitude(s) for diagram number 3789
    FFV1_0( w_fp[516], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3790 OF 15495 ***
    // Wavefunction(s) for diagram number 3790
    FFV1P0_3( w_fp[179], w_fp[259], COUPs[1], 1.0, depCoup, 0., 0., w_fp[517] );
    // Amplitude(s) for diagram number 3790
    VVV1_0( w_fp[451], w_fp[517], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3791 OF 15495 ***
    // Wavefunction(s) for diagram number 3791
    // (none)
    // Amplitude(s) for diagram number 3791
    FFV1_0( w_fp[179], w_fp[439], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];

    // *** DIAGRAM 3792 OF 15495 ***
    // Wavefunction(s) for diagram number 3792
    // (none)
    // Amplitude(s) for diagram number 3792
    FFV1_0( w_fp[181], w_fp[259], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];

    // *** DIAGRAM 3793 OF 15495 ***
    // Wavefunction(s) for diagram number 3793
    // (none)
    // Amplitude(s) for diagram number 3793
    VVV1_0( w_fp[488], w_fp[517], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3794 OF 15495 ***
    // Wavefunction(s) for diagram number 3794
    // (none)
    // Amplitude(s) for diagram number 3794
    FFV1_0( w_fp[179], w_fp[436], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];

    // *** DIAGRAM 3795 OF 15495 ***
    // Wavefunction(s) for diagram number 3795
    // (none)
    // Amplitude(s) for diagram number 3795
    FFV1_0( w_fp[180], w_fp[259], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];

    // *** DIAGRAM 3796 OF 15495 ***
    // Wavefunction(s) for diagram number 3796
    // (none)
    // Amplitude(s) for diagram number 3796
    FFV1_0( w_fp[181], w_fp[436], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3797 OF 15495 ***
    // Wavefunction(s) for diagram number 3797
    // (none)
    // Amplitude(s) for diagram number 3797
    FFV1_0( w_fp[180], w_fp[439], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3798 OF 15495 ***
    // Wavefunction(s) for diagram number 3798
    // (none)
    // Amplitude(s) for diagram number 3798
    FFV1_0( w_fp[179], w_fp[259], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3799 OF 15495 ***
    // Wavefunction(s) for diagram number 3799
    // (none)
    // Amplitude(s) for diagram number 3799
    FFV1_0( w_fp[179], w_fp[453], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];

    // *** DIAGRAM 3800 OF 15495 ***
    // Wavefunction(s) for diagram number 3800
    // (none)
    // Amplitude(s) for diagram number 3800
    FFV1_0( w_fp[516], w_fp[259], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 435 );
    storeWf( wfs, w_cx, nevt, 437 );
    storeWf( wfs, w_cx, nevt, 438 );
    storeWf( wfs, w_cx, nevt, 440 );
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 444 );
    storeWf( wfs, w_cx, nevt, 445 );
    storeWf( wfs, w_cx, nevt, 447 );
    storeWf( wfs, w_cx, nevt, 448 );
    storeWf( wfs, w_cx, nevt, 449 );
    storeWf( wfs, w_cx, nevt, 450 );
    storeWf( wfs, w_cx, nevt, 451 );
    storeWf( wfs, w_cx, nevt, 452 );
    storeWf( wfs, w_cx, nevt, 453 );
    storeWf( wfs, w_cx, nevt, 454 );
    storeWf( wfs, w_cx, nevt, 455 );
    storeWf( wfs, w_cx, nevt, 456 );
    storeWf( wfs, w_cx, nevt, 458 );
    storeWf( wfs, w_cx, nevt, 459 );
    storeWf( wfs, w_cx, nevt, 461 );
    storeWf( wfs, w_cx, nevt, 462 );
    storeWf( wfs, w_cx, nevt, 464 );
    storeWf( wfs, w_cx, nevt, 468 );
    storeWf( wfs, w_cx, nevt, 469 );
    storeWf( wfs, w_cx, nevt, 470 );
    storeWf( wfs, w_cx, nevt, 471 );
    storeWf( wfs, w_cx, nevt, 472 );
    storeWf( wfs, w_cx, nevt, 474 );
    storeWf( wfs, w_cx, nevt, 475 );
    storeWf( wfs, w_cx, nevt, 478 );
    storeWf( wfs, w_cx, nevt, 479 );
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 482 );
    storeWf( wfs, w_cx, nevt, 485 );
    storeWf( wfs, w_cx, nevt, 486 );
    storeWf( wfs, w_cx, nevt, 487 );
    storeWf( wfs, w_cx, nevt, 488 );
    storeWf( wfs, w_cx, nevt, 490 );
    storeWf( wfs, w_cx, nevt, 492 );
    storeWf( wfs, w_cx, nevt, 495 );
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 505 );
    storeWf( wfs, w_cx, nevt, 509 );
    storeWf( wfs, w_cx, nevt, 510 );
    storeWf( wfs, w_cx, nevt, 511 );
    storeWf( wfs, w_cx, nevt, 514 );
    storeWf( wfs, w_cx, nevt, 515 );
    storeWf( wfs, w_cx, nevt, 516 );
    storeWf( wfs, w_cx, nevt, 517 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
