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
  diagramgroup57( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 143 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 230 );
    retrieveWf( wfs, w_cx, nevt, 231 );
    retrieveWf( wfs, w_cx, nevt, 233 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 293 );
    retrieveWf( wfs, w_cx, nevt, 294 );
    retrieveWf( wfs, w_cx, nevt, 295 );
    retrieveWf( wfs, w_cx, nevt, 304 );
    retrieveWf( wfs, w_cx, nevt, 312 );
    retrieveWf( wfs, w_cx, nevt, 313 );
    retrieveWf( wfs, w_cx, nevt, 314 );
    retrieveWf( wfs, w_cx, nevt, 315 );
    retrieveWf( wfs, w_cx, nevt, 316 );
    retrieveWf( wfs, w_cx, nevt, 317 );
    retrieveWf( wfs, w_cx, nevt, 318 );
    retrieveWf( wfs, w_cx, nevt, 319 );
    retrieveWf( wfs, w_cx, nevt, 320 );
    retrieveWf( wfs, w_cx, nevt, 321 );
    retrieveWf( wfs, w_cx, nevt, 322 );
    retrieveWf( wfs, w_cx, nevt, 323 );
    retrieveWf( wfs, w_cx, nevt, 324 );
    retrieveWf( wfs, w_cx, nevt, 326 );
    retrieveWf( wfs, w_cx, nevt, 327 );
    retrieveWf( wfs, w_cx, nevt, 328 );
    retrieveWf( wfs, w_cx, nevt, 329 );
    retrieveWf( wfs, w_cx, nevt, 330 );
    retrieveWf( wfs, w_cx, nevt, 331 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 493 );
    retrieveWf( wfs, w_cx, nevt, 496 );
    retrieveWf( wfs, w_cx, nevt, 498 );
    retrieveWf( wfs, w_cx, nevt, 499 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 504 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 548 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 598 );
    retrieveWf( wfs, w_cx, nevt, 599 );
    retrieveWf( wfs, w_cx, nevt, 602 );
    retrieveWf( wfs, w_cx, nevt, 603 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 696 );
    retrieveWf( wfs, w_cx, nevt, 698 );
    retrieveWf( wfs, w_cx, nevt, 703 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 707 );
    retrieveWf( wfs, w_cx, nevt, 708 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 711 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
    retrieveWf( wfs, w_cx, nevt, 734 );
    retrieveWf( wfs, w_cx, nevt, 735 );
    retrieveWf( wfs, w_cx, nevt, 743 );
    retrieveWf( wfs, w_cx, nevt, 747 );
    retrieveWf( wfs, w_cx, nevt, 748 );
    retrieveWf( wfs, w_cx, nevt, 749 );
#endif
#endif

    // *** DIAGRAM 11201 OF 15495 ***
    // Wavefunction(s) for diagram number 11201
    // (none)
    // Amplitude(s) for diagram number 11201
    FFV1_0( w_fp[3], w_fp[225], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[225], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[225], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11202 OF 15495 ***
    // Wavefunction(s) for diagram number 11202
    // (none)
    // Amplitude(s) for diagram number 11202
    VVV1_0( w_fp[0], w_fp[315], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[316], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[317], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];

    // *** DIAGRAM 11203 OF 15495 ***
    // Wavefunction(s) for diagram number 11203
    // (none)
    // Amplitude(s) for diagram number 11203
    VVV1_0( w_fp[662], w_fp[230], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11204 OF 15495 ***
    // Wavefunction(s) for diagram number 11204
    // (none)
    // Amplitude(s) for diagram number 11204
    FFV1_0( w_fp[196], w_fp[226], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 11205 OF 15495 ***
    // Wavefunction(s) for diagram number 11205
    // (none)
    // Amplitude(s) for diagram number 11205
    FFV1_0( w_fp[216], w_fp[215], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 11206 OF 15495 ***
    // Wavefunction(s) for diagram number 11206
    // (none)
    // Amplitude(s) for diagram number 11206
    FFV1_0( w_fp[326], w_fp[542], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11207 OF 15495 ***
    // Wavefunction(s) for diagram number 11207
    // (none)
    // Amplitude(s) for diagram number 11207
    FFV1_0( w_fp[196], w_fp[542], w_fp[313], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];

    // *** DIAGRAM 11208 OF 15495 ***
    // Wavefunction(s) for diagram number 11208
    // (none)
    // Amplitude(s) for diagram number 11208
    FFV1_0( w_fp[216], w_fp[542], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11209 OF 15495 ***
    // Wavefunction(s) for diagram number 11209
    // (none)
    // Amplitude(s) for diagram number 11209
    FFV1_0( w_fp[536], w_fp[504], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11210 OF 15495 ***
    // Wavefunction(s) for diagram number 11210
    // (none)
    // Amplitude(s) for diagram number 11210
    FFV1_0( w_fp[536], w_fp[215], w_fp[313], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 11211 OF 15495 ***
    // Wavefunction(s) for diagram number 11211
    // (none)
    // Amplitude(s) for diagram number 11211
    FFV1_0( w_fp[536], w_fp[226], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11212 OF 15495 ***
    // Wavefunction(s) for diagram number 11212
    // (none)
    // Amplitude(s) for diagram number 11212
    FFV1_0( w_fp[216], w_fp[504], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11213 OF 15495 ***
    // Wavefunction(s) for diagram number 11213
    // (none)
    // Amplitude(s) for diagram number 11213
    FFV1_0( w_fp[326], w_fp[226], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11214 OF 15495 ***
    // Wavefunction(s) for diagram number 11214
    // (none)
    // Amplitude(s) for diagram number 11214
    VVV1_0( w_fp[0], w_fp[313], w_fp[230], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11215 OF 15495 ***
    // Wavefunction(s) for diagram number 11215
    // (none)
    // Amplitude(s) for diagram number 11215
    FFV1_0( w_fp[196], w_fp[215], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[215], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[215], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11216 OF 15495 ***
    // Wavefunction(s) for diagram number 11216
    // (none)
    // Amplitude(s) for diagram number 11216
    VVV1_0( w_fp[662], w_fp[231], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11217 OF 15495 ***
    // Wavefunction(s) for diagram number 11217
    // (none)
    // Amplitude(s) for diagram number 11217
    FFV1_0( w_fp[168], w_fp[225], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 11218 OF 15495 ***
    // Wavefunction(s) for diagram number 11218
    // (none)
    // Amplitude(s) for diagram number 11218
    FFV1_0( w_fp[218], w_fp[215], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 11219 OF 15495 ***
    // Wavefunction(s) for diagram number 11219
    // (none)
    // Amplitude(s) for diagram number 11219
    FFV1_0( w_fp[327], w_fp[542], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11220 OF 15495 ***
    // Wavefunction(s) for diagram number 11220
    // (none)
    // Amplitude(s) for diagram number 11220
    FFV1_0( w_fp[168], w_fp[542], w_fp[312], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];

    // *** DIAGRAM 11221 OF 15495 ***
    // Wavefunction(s) for diagram number 11221
    // (none)
    // Amplitude(s) for diagram number 11221
    FFV1_0( w_fp[218], w_fp[542], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11222 OF 15495 ***
    // Wavefunction(s) for diagram number 11222
    // (none)
    // Amplitude(s) for diagram number 11222
    FFV1_0( w_fp[527], w_fp[504], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11223 OF 15495 ***
    // Wavefunction(s) for diagram number 11223
    // (none)
    // Amplitude(s) for diagram number 11223
    FFV1_0( w_fp[527], w_fp[215], w_fp[312], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 11224 OF 15495 ***
    // Wavefunction(s) for diagram number 11224
    // (none)
    // Amplitude(s) for diagram number 11224
    FFV1_0( w_fp[527], w_fp[225], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11225 OF 15495 ***
    // Wavefunction(s) for diagram number 11225
    // (none)
    // Amplitude(s) for diagram number 11225
    FFV1_0( w_fp[218], w_fp[504], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11226 OF 15495 ***
    // Wavefunction(s) for diagram number 11226
    // (none)
    // Amplitude(s) for diagram number 11226
    FFV1_0( w_fp[327], w_fp[225], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11227 OF 15495 ***
    // Wavefunction(s) for diagram number 11227
    // (none)
    // Amplitude(s) for diagram number 11227
    VVV1_0( w_fp[0], w_fp[312], w_fp[231], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11228 OF 15495 ***
    // Wavefunction(s) for diagram number 11228
    // (none)
    // Amplitude(s) for diagram number 11228
    FFV1_0( w_fp[168], w_fp[215], w_fp[548], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[500], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11229 OF 15495 ***
    // Wavefunction(s) for diagram number 11229
    // (none)
    // Amplitude(s) for diagram number 11229
    VVV1_0( w_fp[662], w_fp[228], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11230 OF 15495 ***
    // Wavefunction(s) for diagram number 11230
    // (none)
    // Amplitude(s) for diagram number 11230
    FFV1_0( w_fp[3], w_fp[233], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11231 OF 15495 ***
    // Wavefunction(s) for diagram number 11231
    // (none)
    // Amplitude(s) for diagram number 11231
    FFV1_0( w_fp[221], w_fp[215], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11232 OF 15495 ***
    // Wavefunction(s) for diagram number 11232
    // (none)
    // Amplitude(s) for diagram number 11232
    FFV1_0( w_fp[324], w_fp[542], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];

    // *** DIAGRAM 11233 OF 15495 ***
    // Wavefunction(s) for diagram number 11233
    // (none)
    // Amplitude(s) for diagram number 11233
    FFV1_0( w_fp[3], w_fp[542], w_fp[329], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11234 OF 15495 ***
    // Wavefunction(s) for diagram number 11234
    // (none)
    // Amplitude(s) for diagram number 11234
    FFV1_0( w_fp[221], w_fp[542], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];

    // *** DIAGRAM 11235 OF 15495 ***
    // Wavefunction(s) for diagram number 11235
    // (none)
    // Amplitude(s) for diagram number 11235
    FFV1_0( w_fp[3], w_fp[504], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11236 OF 15495 ***
    // Wavefunction(s) for diagram number 11236
    // (none)
    // Amplitude(s) for diagram number 11236
    FFV1_0( w_fp[324], w_fp[215], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11237 OF 15495 ***
    // Wavefunction(s) for diagram number 11237
    // (none)
    // Amplitude(s) for diagram number 11237
    VVV1_0( w_fp[514], w_fp[304], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

    // *** DIAGRAM 11238 OF 15495 ***
    // Wavefunction(s) for diagram number 11238
    // (none)
    // Amplitude(s) for diagram number 11238
    FFV1_0( w_fp[221], w_fp[504], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 11239 OF 15495 ***
    // Wavefunction(s) for diagram number 11239
    // (none)
    // Amplitude(s) for diagram number 11239
    FFV1_0( w_fp[324], w_fp[233], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];

    // *** DIAGRAM 11240 OF 15495 ***
    // Wavefunction(s) for diagram number 11240
    // (none)
    // Amplitude(s) for diagram number 11240
    VVV1_0( w_fp[0], w_fp[329], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 11241 OF 15495 ***
    // Wavefunction(s) for diagram number 11241
    // (none)
    // Amplitude(s) for diagram number 11241
    FFV1_0( w_fp[3], w_fp[215], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[267], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[268], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11242 OF 15495 ***
    // Wavefunction(s) for diagram number 11242
    // (none)
    // Amplitude(s) for diagram number 11242
    VVVV1_0( w_fp[662], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11243 OF 15495 ***
    // Wavefunction(s) for diagram number 11243
    // (none)
    // Amplitude(s) for diagram number 11243
    VVV1_0( w_fp[239], w_fp[7], w_fp[710], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11244 OF 15495 ***
    // Wavefunction(s) for diagram number 11244
    // (none)
    // Amplitude(s) for diagram number 11244
    VVV1_0( w_fp[239], w_fp[5], w_fp[711], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11245 OF 15495 ***
    // Wavefunction(s) for diagram number 11245
    FFV1_1( w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[504] );
    // Amplitude(s) for diagram number 11245
    FFV1_0( w_fp[216], w_fp[504], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 11246 OF 15495 ***
    // Wavefunction(s) for diagram number 11246
    // (none)
    // Amplitude(s) for diagram number 11246
    FFV1_0( w_fp[216], w_fp[2], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11247 OF 15495 ***
    // Wavefunction(s) for diagram number 11247
    // (none)
    // Amplitude(s) for diagram number 11247
    FFV1_0( w_fp[199], w_fp[504], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];

    // *** DIAGRAM 11248 OF 15495 ***
    // Wavefunction(s) for diagram number 11248
    // (none)
    // Amplitude(s) for diagram number 11248
    FFV1_0( w_fp[199], w_fp[2], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11249 OF 15495 ***
    // Wavefunction(s) for diagram number 11249
    // (none)
    // Amplitude(s) for diagram number 11249
    FFV1_0( w_fp[586], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11250 OF 15495 ***
    // Wavefunction(s) for diagram number 11250
    // (none)
    // Amplitude(s) for diagram number 11250
    FFV1_0( w_fp[45], w_fp[496], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11251 OF 15495 ***
    // Wavefunction(s) for diagram number 11251
    // (none)
    // Amplitude(s) for diagram number 11251
    VVV1_0( w_fp[313], w_fp[7], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11252 OF 15495 ***
    // Wavefunction(s) for diagram number 11252
    // (none)
    // Amplitude(s) for diagram number 11252
    FFV1_0( w_fp[45], w_fp[2], w_fp[313], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];

    // *** DIAGRAM 11253 OF 15495 ***
    // Wavefunction(s) for diagram number 11253
    // (none)
    // Amplitude(s) for diagram number 11253
    VVV1_0( w_fp[314], w_fp[5], w_fp[499], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11254 OF 15495 ***
    // Wavefunction(s) for diagram number 11254
    // (none)
    // Amplitude(s) for diagram number 11254
    FFV1_0( w_fp[586], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];

    // *** DIAGRAM 11255 OF 15495 ***
    // Wavefunction(s) for diagram number 11255
    // (none)
    // Amplitude(s) for diagram number 11255
    FFV1_0( w_fp[536], w_fp[2], w_fp[321], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[322], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[536], w_fp[2], w_fp[323], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11256 OF 15495 ***
    // Wavefunction(s) for diagram number 11256
    FFV1_1( w_fp[496], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[503] );
    // Amplitude(s) for diagram number 11256
    FFV1_0( w_fp[216], w_fp[503], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11257 OF 15495 ***
    // Wavefunction(s) for diagram number 11257
    // (none)
    // Amplitude(s) for diagram number 11257
    FFV1_0( w_fp[588], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11258 OF 15495 ***
    // Wavefunction(s) for diagram number 11258
    // (none)
    // Amplitude(s) for diagram number 11258
    FFV1_0( w_fp[199], w_fp[503], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11259 OF 15495 ***
    // Wavefunction(s) for diagram number 11259
    // (none)
    // Amplitude(s) for diagram number 11259
    FFV1_0( w_fp[532], w_fp[496], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11260 OF 15495 ***
    // Wavefunction(s) for diagram number 11260
    // (none)
    // Amplitude(s) for diagram number 11260
    VVVV1_0( w_fp[0], w_fp[313], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[313], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[313], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];

    // *** DIAGRAM 11261 OF 15495 ***
    // Wavefunction(s) for diagram number 11261
    // (none)
    // Amplitude(s) for diagram number 11261
    VVV1_0( w_fp[239], w_fp[7], w_fp[274], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 11262 OF 15495 ***
    // Wavefunction(s) for diagram number 11262
    // (none)
    // Amplitude(s) for diagram number 11262
    VVV1_0( w_fp[313], w_fp[7], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 11263 OF 15495 ***
    // Wavefunction(s) for diagram number 11263
    // (none)
    // Amplitude(s) for diagram number 11263
    FFV1_0( w_fp[199], w_fp[2], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11264 OF 15495 ***
    // Wavefunction(s) for diagram number 11264
    // (none)
    // Amplitude(s) for diagram number 11264
    FFV1_0( w_fp[532], w_fp[2], w_fp[313], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[205] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];

    // *** DIAGRAM 11265 OF 15495 ***
    // Wavefunction(s) for diagram number 11265
    // (none)
    // Amplitude(s) for diagram number 11265
    VVVV1_0( w_fp[0], w_fp[314], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[314], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[314], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 11266 OF 15495 ***
    // Wavefunction(s) for diagram number 11266
    // (none)
    // Amplitude(s) for diagram number 11266
    VVV1_0( w_fp[239], w_fp[5], w_fp[498], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11267 OF 15495 ***
    // Wavefunction(s) for diagram number 11267
    // (none)
    // Amplitude(s) for diagram number 11267
    VVV1_0( w_fp[314], w_fp[5], w_fp[272], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 11268 OF 15495 ***
    // Wavefunction(s) for diagram number 11268
    // (none)
    // Amplitude(s) for diagram number 11268
    FFV1_0( w_fp[216], w_fp[2], w_fp[498], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11269 OF 15495 ***
    // Wavefunction(s) for diagram number 11269
    // (none)
    // Amplitude(s) for diagram number 11269
    FFV1_0( w_fp[588], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 11270 OF 15495 ***
    // Wavefunction(s) for diagram number 11270
    // (none)
    // Amplitude(s) for diagram number 11270
    VVV1_0( w_fp[670], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVV1_0( w_fp[703], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVV1_0( w_fp[709], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 11271 OF 15495 ***
    // Wavefunction(s) for diagram number 11271
    // (none)
    // Amplitude(s) for diagram number 11271
    FFV1_0( w_fp[199], w_fp[2], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11272 OF 15495 ***
    // Wavefunction(s) for diagram number 11272
    // (none)
    // Amplitude(s) for diagram number 11272
    VVV1_0( w_fp[748], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVV1_0( w_fp[708], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    VVV1_0( w_fp[712], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 11273 OF 15495 ***
    // Wavefunction(s) for diagram number 11273
    // (none)
    // Amplitude(s) for diagram number 11273
    FFV1_0( w_fp[216], w_fp[2], w_fp[748], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11274 OF 15495 ***
    // Wavefunction(s) for diagram number 11274
    // (none)
    // Amplitude(s) for diagram number 11274
    VVV1_0( w_fp[0], w_fp[321], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[322], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[323], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 11275 OF 15495 ***
    // Wavefunction(s) for diagram number 11275
    // (none)
    // Amplitude(s) for diagram number 11275
    VVV1_0( w_fp[662], w_fp[239], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 11276 OF 15495 ***
    // Wavefunction(s) for diagram number 11276
    // (none)
    // Amplitude(s) for diagram number 11276
    FFV1_0( w_fp[196], w_fp[122], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11277 OF 15495 ***
    // Wavefunction(s) for diagram number 11277
    // (none)
    // Amplitude(s) for diagram number 11277
    FFV1_0( w_fp[123], w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11278 OF 15495 ***
    // Wavefunction(s) for diagram number 11278
    // (none)
    // Amplitude(s) for diagram number 11278
    FFV1_0( w_fp[536], w_fp[496], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];

    // *** DIAGRAM 11279 OF 15495 ***
    // Wavefunction(s) for diagram number 11279
    // (none)
    // Amplitude(s) for diagram number 11279
    FFV1_0( w_fp[536], w_fp[2], w_fp[331], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11280 OF 15495 ***
    // Wavefunction(s) for diagram number 11280
    // (none)
    // Amplitude(s) for diagram number 11280
    FFV1_0( w_fp[536], w_fp[122], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[466] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 11281 OF 15495 ***
    // Wavefunction(s) for diagram number 11281
    // (none)
    // Amplitude(s) for diagram number 11281
    FFV1_0( w_fp[196], w_fp[496], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11282 OF 15495 ***
    // Wavefunction(s) for diagram number 11282
    // (none)
    // Amplitude(s) for diagram number 11282
    FFV1_0( w_fp[326], w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11283 OF 15495 ***
    // Wavefunction(s) for diagram number 11283
    // (none)
    // Amplitude(s) for diagram number 11283
    VVV1_0( w_fp[522], w_fp[304], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];

    // *** DIAGRAM 11284 OF 15495 ***
    // Wavefunction(s) for diagram number 11284
    // (none)
    // Amplitude(s) for diagram number 11284
    FFV1_0( w_fp[123], w_fp[496], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 11285 OF 15495 ***
    // Wavefunction(s) for diagram number 11285
    // (none)
    // Amplitude(s) for diagram number 11285
    FFV1_0( w_fp[326], w_fp[122], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];

    // *** DIAGRAM 11286 OF 15495 ***
    // Wavefunction(s) for diagram number 11286
    // (none)
    // Amplitude(s) for diagram number 11286
    VVV1_0( w_fp[0], w_fp[331], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];

    // *** DIAGRAM 11287 OF 15495 ***
    // Wavefunction(s) for diagram number 11287
    // (none)
    // Amplitude(s) for diagram number 11287
    FFV1_0( w_fp[196], w_fp[2], w_fp[493], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[714], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[196], w_fp[2], w_fp[713], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 11288 OF 15495 ***
    // Wavefunction(s) for diagram number 11288
    // (none)
    // Amplitude(s) for diagram number 11288
    VVVV1_0( w_fp[662], w_fp[154], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[154], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[154], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 11289 OF 15495 ***
    // Wavefunction(s) for diagram number 11289
    // (none)
    // Amplitude(s) for diagram number 11289
    VVV1_0( w_fp[154], w_fp[7], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 11290 OF 15495 ***
    // Wavefunction(s) for diagram number 11290
    // (none)
    // Amplitude(s) for diagram number 11290
    VVV1_0( w_fp[154], w_fp[4], w_fp[711], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 11291 OF 15495 ***
    // Wavefunction(s) for diagram number 11291
    // (none)
    // Amplitude(s) for diagram number 11291
    FFV1_0( w_fp[218], w_fp[504], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];

    // *** DIAGRAM 11292 OF 15495 ***
    // Wavefunction(s) for diagram number 11292
    // (none)
    // Amplitude(s) for diagram number 11292
    FFV1_0( w_fp[218], w_fp[2], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11293 OF 15495 ***
    // Wavefunction(s) for diagram number 11293
    // (none)
    // Amplitude(s) for diagram number 11293
    FFV1_0( w_fp[171], w_fp[504], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];

    // *** DIAGRAM 11294 OF 15495 ***
    // Wavefunction(s) for diagram number 11294
    // (none)
    // Amplitude(s) for diagram number 11294
    FFV1_0( w_fp[171], w_fp[2], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11295 OF 15495 ***
    // Wavefunction(s) for diagram number 11295
    // (none)
    // Amplitude(s) for diagram number 11295
    FFV1_0( w_fp[509], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11296 OF 15495 ***
    // Wavefunction(s) for diagram number 11296
    // (none)
    // Amplitude(s) for diagram number 11296
    FFV1_0( w_fp[543], w_fp[496], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11297 OF 15495 ***
    // Wavefunction(s) for diagram number 11297
    // (none)
    // Amplitude(s) for diagram number 11297
    VVV1_0( w_fp[312], w_fp[7], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11298 OF 15495 ***
    // Wavefunction(s) for diagram number 11298
    // (none)
    // Amplitude(s) for diagram number 11298
    FFV1_0( w_fp[543], w_fp[2], w_fp[312], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];

    // *** DIAGRAM 11299 OF 15495 ***
    // Wavefunction(s) for diagram number 11299
    // (none)
    // Amplitude(s) for diagram number 11299
    VVV1_0( w_fp[314], w_fp[4], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11300 OF 15495 ***
    // Wavefunction(s) for diagram number 11300
    // (none)
    // Amplitude(s) for diagram number 11300
    FFV1_0( w_fp[509], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];

    // *** DIAGRAM 11301 OF 15495 ***
    // Wavefunction(s) for diagram number 11301
    // (none)
    // Amplitude(s) for diagram number 11301
    FFV1_0( w_fp[527], w_fp[2], w_fp[318], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[319], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[320], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11302 OF 15495 ***
    // Wavefunction(s) for diagram number 11302
    // (none)
    // Amplitude(s) for diagram number 11302
    FFV1_0( w_fp[218], w_fp[503], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11303 OF 15495 ***
    // Wavefunction(s) for diagram number 11303
    // (none)
    // Amplitude(s) for diagram number 11303
    FFV1_0( w_fp[481], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11304 OF 15495 ***
    // Wavefunction(s) for diagram number 11304
    // (none)
    // Amplitude(s) for diagram number 11304
    FFV1_0( w_fp[171], w_fp[503], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11305 OF 15495 ***
    // Wavefunction(s) for diagram number 11305
    // (none)
    // Amplitude(s) for diagram number 11305
    FFV1_0( w_fp[450], w_fp[496], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11306 OF 15495 ***
    // Wavefunction(s) for diagram number 11306
    // (none)
    // Amplitude(s) for diagram number 11306
    VVVV1_0( w_fp[0], w_fp[312], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[312], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[312], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];

    // *** DIAGRAM 11307 OF 15495 ***
    // Wavefunction(s) for diagram number 11307
    // (none)
    // Amplitude(s) for diagram number 11307
    VVV1_0( w_fp[154], w_fp[7], w_fp[743], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 11308 OF 15495 ***
    // Wavefunction(s) for diagram number 11308
    // (none)
    // Amplitude(s) for diagram number 11308
    VVV1_0( w_fp[312], w_fp[7], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 11309 OF 15495 ***
    // Wavefunction(s) for diagram number 11309
    // (none)
    // Amplitude(s) for diagram number 11309
    FFV1_0( w_fp[171], w_fp[2], w_fp[743], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11310 OF 15495 ***
    // Wavefunction(s) for diagram number 11310
    // (none)
    // Amplitude(s) for diagram number 11310
    FFV1_0( w_fp[450], w_fp[2], w_fp[312], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[199] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];

    // *** DIAGRAM 11311 OF 15495 ***
    // Wavefunction(s) for diagram number 11311
    // (none)
    // Amplitude(s) for diagram number 11311
    VVVV1_0( w_fp[0], w_fp[314], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[314], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[314], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 11312 OF 15495 ***
    // Wavefunction(s) for diagram number 11312
    // (none)
    // Amplitude(s) for diagram number 11312
    VVV1_0( w_fp[154], w_fp[4], w_fp[498], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 11313 OF 15495 ***
    // Wavefunction(s) for diagram number 11313
    // (none)
    // Amplitude(s) for diagram number 11313
    VVV1_0( w_fp[314], w_fp[4], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 11314 OF 15495 ***
    // Wavefunction(s) for diagram number 11314
    // (none)
    // Amplitude(s) for diagram number 11314
    FFV1_0( w_fp[218], w_fp[2], w_fp[498], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11315 OF 15495 ***
    // Wavefunction(s) for diagram number 11315
    // (none)
    // Amplitude(s) for diagram number 11315
    FFV1_0( w_fp[481], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 11316 OF 15495 ***
    // Wavefunction(s) for diagram number 11316
    // (none)
    // Amplitude(s) for diagram number 11316
    VVV1_0( w_fp[548], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVV1_0( w_fp[500], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVV1_0( w_fp[706], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 11317 OF 15495 ***
    // Wavefunction(s) for diagram number 11317
    // (none)
    // Amplitude(s) for diagram number 11317
    FFV1_0( w_fp[171], w_fp[2], w_fp[548], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[500], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11318 OF 15495 ***
    // Wavefunction(s) for diagram number 11318
    // (none)
    // Amplitude(s) for diagram number 11318
    VVV1_0( w_fp[748], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVV1_0( w_fp[708], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    VVV1_0( w_fp[712], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 11319 OF 15495 ***
    // Wavefunction(s) for diagram number 11319
    // (none)
    // Amplitude(s) for diagram number 11319
    FFV1_0( w_fp[218], w_fp[2], w_fp[748], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11320 OF 15495 ***
    // Wavefunction(s) for diagram number 11320
    // (none)
    // Amplitude(s) for diagram number 11320
    VVV1_0( w_fp[0], w_fp[318], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[319], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[320], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 11321 OF 15495 ***
    // Wavefunction(s) for diagram number 11321
    // (none)
    // Amplitude(s) for diagram number 11321
    VVV1_0( w_fp[662], w_fp[154], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 11322 OF 15495 ***
    // Wavefunction(s) for diagram number 11322
    // (none)
    // Amplitude(s) for diagram number 11322
    FFV1_0( w_fp[168], w_fp[98], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11323 OF 15495 ***
    // Wavefunction(s) for diagram number 11323
    // (none)
    // Amplitude(s) for diagram number 11323
    FFV1_0( w_fp[99], w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11324 OF 15495 ***
    // Wavefunction(s) for diagram number 11324
    // (none)
    // Amplitude(s) for diagram number 11324
    FFV1_0( w_fp[527], w_fp[496], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];

    // *** DIAGRAM 11325 OF 15495 ***
    // Wavefunction(s) for diagram number 11325
    // (none)
    // Amplitude(s) for diagram number 11325
    FFV1_0( w_fp[527], w_fp[2], w_fp[330], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11326 OF 15495 ***
    // Wavefunction(s) for diagram number 11326
    // (none)
    // Amplitude(s) for diagram number 11326
    FFV1_0( w_fp[527], w_fp[98], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 11327 OF 15495 ***
    // Wavefunction(s) for diagram number 11327
    // (none)
    // Amplitude(s) for diagram number 11327
    FFV1_0( w_fp[168], w_fp[496], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11328 OF 15495 ***
    // Wavefunction(s) for diagram number 11328
    // (none)
    // Amplitude(s) for diagram number 11328
    FFV1_0( w_fp[327], w_fp[2], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11329 OF 15495 ***
    // Wavefunction(s) for diagram number 11329
    // (none)
    // Amplitude(s) for diagram number 11329
    VVV1_0( w_fp[518], w_fp[304], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 11330 OF 15495 ***
    // Wavefunction(s) for diagram number 11330
    // (none)
    // Amplitude(s) for diagram number 11330
    FFV1_0( w_fp[99], w_fp[496], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[193] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];

    // *** DIAGRAM 11331 OF 15495 ***
    // Wavefunction(s) for diagram number 11331
    // (none)
    // Amplitude(s) for diagram number 11331
    FFV1_0( w_fp[327], w_fp[98], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 11332 OF 15495 ***
    // Wavefunction(s) for diagram number 11332
    // (none)
    // Amplitude(s) for diagram number 11332
    VVV1_0( w_fp[0], w_fp[330], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];

    // *** DIAGRAM 11333 OF 15495 ***
    // Wavefunction(s) for diagram number 11333
    // (none)
    // Amplitude(s) for diagram number 11333
    FFV1_0( w_fp[168], w_fp[2], w_fp[295], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[73] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[346] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[294], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[168], w_fp[2], w_fp[293], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[346] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 11334 OF 15495 ***
    // Wavefunction(s) for diagram number 11334
    // (none)
    // Amplitude(s) for diagram number 11334
    VVVV1_0( w_fp[662], w_fp[144], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[144], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[144], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 11335 OF 15495 ***
    // Wavefunction(s) for diagram number 11335
    // (none)
    // Amplitude(s) for diagram number 11335
    VVV1_0( w_fp[144], w_fp[5], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];

    // *** DIAGRAM 11336 OF 15495 ***
    // Wavefunction(s) for diagram number 11336
    // (none)
    // Amplitude(s) for diagram number 11336
    VVV1_0( w_fp[144], w_fp[4], w_fp[710], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 11337 OF 15495 ***
    // Wavefunction(s) for diagram number 11337
    // (none)
    // Amplitude(s) for diagram number 11337
    FFV1_0( w_fp[204], w_fp[504], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 11338 OF 15495 ***
    // Wavefunction(s) for diagram number 11338
    // (none)
    // Amplitude(s) for diagram number 11338
    FFV1_0( w_fp[204], w_fp[2], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11339 OF 15495 ***
    // Wavefunction(s) for diagram number 11339
    // (none)
    // Amplitude(s) for diagram number 11339
    FFV1_0( w_fp[180], w_fp[504], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];

    // *** DIAGRAM 11340 OF 15495 ***
    // Wavefunction(s) for diagram number 11340
    // (none)
    // Amplitude(s) for diagram number 11340
    FFV1_0( w_fp[180], w_fp[2], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11341 OF 15495 ***
    // Wavefunction(s) for diagram number 11341
    // (none)
    // Amplitude(s) for diagram number 11341
    FFV1_0( w_fp[435], w_fp[496], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11342 OF 15495 ***
    // Wavefunction(s) for diagram number 11342
    // (none)
    // Amplitude(s) for diagram number 11342
    FFV1_0( w_fp[438], w_fp[496], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11343 OF 15495 ***
    // Wavefunction(s) for diagram number 11343
    // (none)
    // Amplitude(s) for diagram number 11343
    VVV1_0( w_fp[312], w_fp[5], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11344 OF 15495 ***
    // Wavefunction(s) for diagram number 11344
    // (none)
    // Amplitude(s) for diagram number 11344
    FFV1_0( w_fp[438], w_fp[2], w_fp[312], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];

    // *** DIAGRAM 11345 OF 15495 ***
    // Wavefunction(s) for diagram number 11345
    // (none)
    // Amplitude(s) for diagram number 11345
    VVV1_0( w_fp[313], w_fp[4], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11346 OF 15495 ***
    // Wavefunction(s) for diagram number 11346
    // (none)
    // Amplitude(s) for diagram number 11346
    FFV1_0( w_fp[435], w_fp[2], w_fp[313], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];

    // *** DIAGRAM 11347 OF 15495 ***
    // Wavefunction(s) for diagram number 11347
    // (none)
    // Amplitude(s) for diagram number 11347
    FFV1_0( w_fp[596], w_fp[2], w_fp[315], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[316], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[317], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11348 OF 15495 ***
    // Wavefunction(s) for diagram number 11348
    // (none)
    // Amplitude(s) for diagram number 11348
    FFV1_0( w_fp[204], w_fp[503], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11349 OF 15495 ***
    // Wavefunction(s) for diagram number 11349
    // (none)
    // Amplitude(s) for diagram number 11349
    FFV1_0( w_fp[584], w_fp[496], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11350 OF 15495 ***
    // Wavefunction(s) for diagram number 11350
    // (none)
    // Amplitude(s) for diagram number 11350
    FFV1_0( w_fp[180], w_fp[503], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11351 OF 15495 ***
    // Wavefunction(s) for diagram number 11351
    // (none)
    // Amplitude(s) for diagram number 11351
    FFV1_0( w_fp[557], w_fp[496], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11352 OF 15495 ***
    // Wavefunction(s) for diagram number 11352
    // (none)
    // Amplitude(s) for diagram number 11352
    VVVV1_0( w_fp[0], w_fp[312], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[312], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[312], w_fp[144], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];

    // *** DIAGRAM 11353 OF 15495 ***
    // Wavefunction(s) for diagram number 11353
    // (none)
    // Amplitude(s) for diagram number 11353
    VVV1_0( w_fp[144], w_fp[5], w_fp[743], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];

    // *** DIAGRAM 11354 OF 15495 ***
    // Wavefunction(s) for diagram number 11354
    // (none)
    // Amplitude(s) for diagram number 11354
    VVV1_0( w_fp[312], w_fp[5], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];

    // *** DIAGRAM 11355 OF 15495 ***
    // Wavefunction(s) for diagram number 11355
    // (none)
    // Amplitude(s) for diagram number 11355
    FFV1_0( w_fp[180], w_fp[2], w_fp[743], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11356 OF 15495 ***
    // Wavefunction(s) for diagram number 11356
    // (none)
    // Amplitude(s) for diagram number 11356
    FFV1_0( w_fp[557], w_fp[2], w_fp[312], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[198] += amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];

    // *** DIAGRAM 11357 OF 15495 ***
    // Wavefunction(s) for diagram number 11357
    // (none)
    // Amplitude(s) for diagram number 11357
    VVVV1_0( w_fp[0], w_fp[313], w_fp[144], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[313], w_fp[144], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[313], w_fp[144], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];

    // *** DIAGRAM 11358 OF 15495 ***
    // Wavefunction(s) for diagram number 11358
    // (none)
    // Amplitude(s) for diagram number 11358
    VVV1_0( w_fp[144], w_fp[4], w_fp[274], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];

    // *** DIAGRAM 11359 OF 15495 ***
    // Wavefunction(s) for diagram number 11359
    // (none)
    // Amplitude(s) for diagram number 11359
    VVV1_0( w_fp[313], w_fp[4], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];

    // *** DIAGRAM 11360 OF 15495 ***
    // Wavefunction(s) for diagram number 11360
    // (none)
    // Amplitude(s) for diagram number 11360
    FFV1_0( w_fp[204], w_fp[2], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11361 OF 15495 ***
    // Wavefunction(s) for diagram number 11361
    // (none)
    // Amplitude(s) for diagram number 11361
    FFV1_0( w_fp[584], w_fp[2], w_fp[313], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[204] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];

    // *** DIAGRAM 11362 OF 15495 ***
    // Wavefunction(s) for diagram number 11362
    // (none)
    // Amplitude(s) for diagram number 11362
    VVV1_0( w_fp[548], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[276] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[372] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    VVV1_0( w_fp[500], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    VVV1_0( w_fp[706], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];

    // *** DIAGRAM 11363 OF 15495 ***
    // Wavefunction(s) for diagram number 11363
    // (none)
    // Amplitude(s) for diagram number 11363
    FFV1_0( w_fp[180], w_fp[2], w_fp[548], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[180], w_fp[2], w_fp[500], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[180], w_fp[2], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11364 OF 15495 ***
    // Wavefunction(s) for diagram number 11364
    // (none)
    // Amplitude(s) for diagram number 11364
    VVV1_0( w_fp[670], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    VVV1_0( w_fp[703], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    VVV1_0( w_fp[709], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[276] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];

    // *** DIAGRAM 11365 OF 15495 ***
    // Wavefunction(s) for diagram number 11365
    // (none)
    // Amplitude(s) for diagram number 11365
    FFV1_0( w_fp[204], w_fp[2], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[204], w_fp[2], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[204], w_fp[2], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11366 OF 15495 ***
    // Wavefunction(s) for diagram number 11366
    // (none)
    // Amplitude(s) for diagram number 11366
    VVV1_0( w_fp[0], w_fp[315], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[316], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[317], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[278] += amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];

    // *** DIAGRAM 11367 OF 15495 ***
    // Wavefunction(s) for diagram number 11367
    // (none)
    // Amplitude(s) for diagram number 11367
    VVV1_0( w_fp[662], w_fp[144], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];

    // *** DIAGRAM 11368 OF 15495 ***
    // Wavefunction(s) for diagram number 11368
    // (none)
    // Amplitude(s) for diagram number 11368
    FFV1_0( w_fp[179], w_fp[148], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11369 OF 15495 ***
    // Wavefunction(s) for diagram number 11369
    // (none)
    // Amplitude(s) for diagram number 11369
    FFV1_0( w_fp[143], w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11370 OF 15495 ***
    // Wavefunction(s) for diagram number 11370
    // (none)
    // Amplitude(s) for diagram number 11370
    FFV1_0( w_fp[596], w_fp[496], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];

    // *** DIAGRAM 11371 OF 15495 ***
    // Wavefunction(s) for diagram number 11371
    // (none)
    // Amplitude(s) for diagram number 11371
    FFV1_0( w_fp[596], w_fp[2], w_fp[329], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11372 OF 15495 ***
    // Wavefunction(s) for diagram number 11372
    // (none)
    // Amplitude(s) for diagram number 11372
    FFV1_0( w_fp[596], w_fp[148], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];

    // *** DIAGRAM 11373 OF 15495 ***
    // Wavefunction(s) for diagram number 11373
    // (none)
    // Amplitude(s) for diagram number 11373
    FFV1_0( w_fp[179], w_fp[496], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11374 OF 15495 ***
    // Wavefunction(s) for diagram number 11374
    // (none)
    // Amplitude(s) for diagram number 11374
    FFV1_0( w_fp[328], w_fp[2], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11375 OF 15495 ***
    // Wavefunction(s) for diagram number 11375
    // (none)
    // Amplitude(s) for diagram number 11375
    VVV1_0( w_fp[514], w_fp[304], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];

    // *** DIAGRAM 11376 OF 15495 ***
    // Wavefunction(s) for diagram number 11376
    // (none)
    // Amplitude(s) for diagram number 11376
    FFV1_0( w_fp[143], w_fp[496], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 11377 OF 15495 ***
    // Wavefunction(s) for diagram number 11377
    // (none)
    // Amplitude(s) for diagram number 11377
    FFV1_0( w_fp[328], w_fp[148], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];

    // *** DIAGRAM 11378 OF 15495 ***
    // Wavefunction(s) for diagram number 11378
    // (none)
    // Amplitude(s) for diagram number 11378
    VVV1_0( w_fp[0], w_fp[329], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];

    // *** DIAGRAM 11379 OF 15495 ***
    // Wavefunction(s) for diagram number 11379
    // (none)
    // Amplitude(s) for diagram number 11379
    FFV1_0( w_fp[179], w_fp[2], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[416] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[267], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[268], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];

    // *** DIAGRAM 11380 OF 15495 ***
    // Wavefunction(s) for diagram number 11380
    // (none)
    // Amplitude(s) for diagram number 11380
    FFV1_0( w_fp[551], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11381 OF 15495 ***
    // Wavefunction(s) for diagram number 11381
    // (none)
    // Amplitude(s) for diagram number 11381
    FFV1_0( w_fp[3], w_fp[148], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];

    // *** DIAGRAM 11382 OF 15495 ***
    // Wavefunction(s) for diagram number 11382
    // (none)
    // Amplitude(s) for diagram number 11382
    FFV1_0( w_fp[221], w_fp[504], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11383 OF 15495 ***
    // Wavefunction(s) for diagram number 11383
    // (none)
    // Amplitude(s) for diagram number 11383
    FFV1_0( w_fp[221], w_fp[2], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 11384 OF 15495 ***
    // Wavefunction(s) for diagram number 11384
    // (none)
    // Amplitude(s) for diagram number 11384
    FFV1_0( w_fp[3], w_fp[504], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];

    // *** DIAGRAM 11385 OF 15495 ***
    // Wavefunction(s) for diagram number 11385
    // (none)
    // Amplitude(s) for diagram number 11385
    FFV1_0( w_fp[551], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[309] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 11386 OF 15495 ***
    // Wavefunction(s) for diagram number 11386
    // (none)
    // Amplitude(s) for diagram number 11386
    FFV1_0( w_fp[599], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11387 OF 15495 ***
    // Wavefunction(s) for diagram number 11387
    // (none)
    // Amplitude(s) for diagram number 11387
    FFV1_0( w_fp[3], w_fp[496], w_fp[598], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];

    // *** DIAGRAM 11388 OF 15495 ***
    // Wavefunction(s) for diagram number 11388
    // (none)
    // Amplitude(s) for diagram number 11388
    FFV1_0( w_fp[324], w_fp[696], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11389 OF 15495 ***
    // Wavefunction(s) for diagram number 11389
    // (none)
    // Amplitude(s) for diagram number 11389
    FFV1_0( w_fp[324], w_fp[2], w_fp[598], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];

    // *** DIAGRAM 11390 OF 15495 ***
    // Wavefunction(s) for diagram number 11390
    // (none)
    // Amplitude(s) for diagram number 11390
    FFV1_0( w_fp[3], w_fp[696], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];

    // *** DIAGRAM 11391 OF 15495 ***
    // Wavefunction(s) for diagram number 11391
    // (none)
    // Amplitude(s) for diagram number 11391
    FFV1_0( w_fp[599], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

    // *** DIAGRAM 11392 OF 15495 ***
    // Wavefunction(s) for diagram number 11392
    // (none)
    // Amplitude(s) for diagram number 11392
    FFV1_0( w_fp[221], w_fp[503], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 11393 OF 15495 ***
    // Wavefunction(s) for diagram number 11393
    // (none)
    // Amplitude(s) for diagram number 11393
    FFV1_0( w_fp[603], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];

    // *** DIAGRAM 11394 OF 15495 ***
    // Wavefunction(s) for diagram number 11394
    // (none)
    // Amplitude(s) for diagram number 11394
    FFV1_0( w_fp[3], w_fp[503], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11395 OF 15495 ***
    // Wavefunction(s) for diagram number 11395
    // (none)
    // Amplitude(s) for diagram number 11395
    FFV1_0( w_fp[3], w_fp[496], w_fp[602], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];

    // *** DIAGRAM 11396 OF 15495 ***
    // Wavefunction(s) for diagram number 11396
    // (none)
    // Amplitude(s) for diagram number 11396
    FFV1_0( w_fp[270], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];

    // *** DIAGRAM 11397 OF 15495 ***
    // Wavefunction(s) for diagram number 11397
    // (none)
    // Amplitude(s) for diagram number 11397
    FFV1_0( w_fp[324], w_fp[698], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];

    // *** DIAGRAM 11398 OF 15495 ***
    // Wavefunction(s) for diagram number 11398
    // (none)
    // Amplitude(s) for diagram number 11398
    FFV1_0( w_fp[270], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11399 OF 15495 ***
    // Wavefunction(s) for diagram number 11399
    // (none)
    // Amplitude(s) for diagram number 11399
    FFV1_0( w_fp[324], w_fp[2], w_fp[602], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];

    // *** DIAGRAM 11400 OF 15495 ***
    // Wavefunction(s) for diagram number 11400
    // (none)
    // Amplitude(s) for diagram number 11400
    FFV1_0( w_fp[3], w_fp[148], w_fp[498], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 503 );
    storeWf( wfs, w_cx, nevt, 504 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
