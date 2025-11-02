// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessChannelIds.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessWavefunctions.h"
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
  diagramgroup33( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                  fptype* jamps,                  // output jamps[ncolor*2*nevt]
                  const int nGoodHel,             // input: number of good helicities
                  const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                  cxtype_sv* jamps,               // output jamps[ncolor*2*neppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 209 );
    retrieveWf( wfs, w_cx, nevt, 210 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 222 );
    retrieveWf( wfs, w_cx, nevt, 223 );
    retrieveWf( wfs, w_cx, nevt, 224 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 440 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 446 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 455 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 458 );
    retrieveWf( wfs, w_cx, nevt, 461 );
    retrieveWf( wfs, w_cx, nevt, 462 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 478 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 511 );
    retrieveWf( wfs, w_cx, nevt, 513 );
    retrieveWf( wfs, w_cx, nevt, 515 );
#endif
#endif

    // *** DIAGRAM 3201 OF 15495 ***
    // Wavefunction(s) for diagram number 3201
    // (none)
    // Amplitude(s) for diagram number 3201
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[169], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];

    // *** DIAGRAM 3202 OF 15495 ***
    // Wavefunction(s) for diagram number 3202
    // (none)
    // Amplitude(s) for diagram number 3202
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[363], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3203 OF 15495 ***
    // Wavefunction(s) for diagram number 3203
    // (none)
    // Amplitude(s) for diagram number 3203
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[104], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3204 OF 15495 ***
    // Wavefunction(s) for diagram number 3204
    // (none)
    // Amplitude(s) for diagram number 3204
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[6], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[6], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[6], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3205 OF 15495 ***
    // Wavefunction(s) for diagram number 3205
    // (none)
    // Amplitude(s) for diagram number 3205
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[462], w_fp[477], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3206 OF 15495 ***
    // Wavefunction(s) for diagram number 3206
    // (none)
    // Amplitude(s) for diagram number 3206
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[462], w_fp[191], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3207 OF 15495 ***
    // Wavefunction(s) for diagram number 3207
    // (none)
    // Amplitude(s) for diagram number 3207
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[477], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3208 OF 15495 ***
    // Wavefunction(s) for diagram number 3208
    // (none)
    // Amplitude(s) for diagram number 3208
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[169], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[429] -= amp_sv[0];
    jamp_sv[471] += amp_sv[0];

    // *** DIAGRAM 3209 OF 15495 ***
    // Wavefunction(s) for diagram number 3209
    // (none)
    // Amplitude(s) for diagram number 3209
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[447], w_fp[209], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3210 OF 15495 ***
    // Wavefunction(s) for diagram number 3210
    // (none)
    // Amplitude(s) for diagram number 3210
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[477], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];

    // *** DIAGRAM 3211 OF 15495 ***
    // Wavefunction(s) for diagram number 3211
    // (none)
    // Amplitude(s) for diagram number 3211
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[191], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];

    // *** DIAGRAM 3212 OF 15495 ***
    // Wavefunction(s) for diagram number 3212
    // (none)
    // Amplitude(s) for diagram number 3212
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[190], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3213 OF 15495 ***
    // Wavefunction(s) for diagram number 3213
    // (none)
    // Amplitude(s) for diagram number 3213
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[210], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3214 OF 15495 ***
    // Wavefunction(s) for diagram number 3214
    // (none)
    // Amplitude(s) for diagram number 3214
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[169], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 3215 OF 15495 ***
    // Wavefunction(s) for diagram number 3215
    // (none)
    // Amplitude(s) for diagram number 3215
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[361], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3216 OF 15495 ***
    // Wavefunction(s) for diagram number 3216
    // (none)
    // Amplitude(s) for diagram number 3216
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[130], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3217 OF 15495 ***
    // Wavefunction(s) for diagram number 3217
    // (none)
    // Amplitude(s) for diagram number 3217
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[84], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[84], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[84], w_fp[488], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3218 OF 15495 ***
    // Wavefunction(s) for diagram number 3218
    // (none)
    // Amplitude(s) for diagram number 3218
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[477], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3219 OF 15495 ***
    // Wavefunction(s) for diagram number 3219
    // (none)
    // Amplitude(s) for diagram number 3219
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[169], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];

    // *** DIAGRAM 3220 OF 15495 ***
    // Wavefunction(s) for diagram number 3220
    // (none)
    // Amplitude(s) for diagram number 3220
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[210], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3221 OF 15495 ***
    // Wavefunction(s) for diagram number 3221
    // (none)
    // Amplitude(s) for diagram number 3221
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[477], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3222 OF 15495 ***
    // Wavefunction(s) for diagram number 3222
    // (none)
    // Amplitude(s) for diagram number 3222
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[458], w_fp[190], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3223 OF 15495 ***
    // Wavefunction(s) for diagram number 3223
    // (none)
    // Amplitude(s) for diagram number 3223
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[477], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 3224 OF 15495 ***
    // Wavefunction(s) for diagram number 3224
    // (none)
    // Amplitude(s) for diagram number 3224
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[190], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];

    // *** DIAGRAM 3225 OF 15495 ***
    // Wavefunction(s) for diagram number 3225
    // (none)
    // Amplitude(s) for diagram number 3225
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[169], w_fp[145], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[169], w_fp[146], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[169], w_fp[147], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 3226 OF 15495 ***
    // Wavefunction(s) for diagram number 3226
    // (none)
    // Amplitude(s) for diagram number 3226
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[145], w_fp[488], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[146], w_fp[488], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[147], w_fp[488], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3227 OF 15495 ***
    // Wavefunction(s) for diagram number 3227
    // (none)
    // Amplitude(s) for diagram number 3227
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[169], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[475], w_fp[169], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[476], w_fp[169], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];

    // *** DIAGRAM 3228 OF 15495 ***
    // Wavefunction(s) for diagram number 3228
    // (none)
    // Amplitude(s) for diagram number 3228
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[511], w_fp[211], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[551] -= amp_sv[0];

    // *** DIAGRAM 3229 OF 15495 ***
    // Wavefunction(s) for diagram number 3229
    // (none)
    // Amplitude(s) for diagram number 3229
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[211], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] -= amp_sv[0];

    // *** DIAGRAM 3230 OF 15495 ***
    // Wavefunction(s) for diagram number 3230
    // (none)
    // Amplitude(s) for diagram number 3230
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[212], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[575] -= amp_sv[0];

    // *** DIAGRAM 3231 OF 15495 ***
    // Wavefunction(s) for diagram number 3231
    // (none)
    // Amplitude(s) for diagram number 3231
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[212], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[569] -= amp_sv[0];

    // *** DIAGRAM 3232 OF 15495 ***
    // Wavefunction(s) for diagram number 3232
    // (none)
    // Amplitude(s) for diagram number 3232
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[213], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 3233 OF 15495 ***
    // Wavefunction(s) for diagram number 3233
    // (none)
    // Amplitude(s) for diagram number 3233
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[511], w_fp[213], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[593] -= amp_sv[0];

    // *** DIAGRAM 3234 OF 15495 ***
    // Wavefunction(s) for diagram number 3234
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[197], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[476] );
    // Amplitude(s) for diagram number 3234
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[440], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[527] -= amp_sv[0];

    // *** DIAGRAM 3235 OF 15495 ***
    // Wavefunction(s) for diagram number 3235
    // (none)
    // Amplitude(s) for diagram number 3235
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] -= amp_sv[0];

    // *** DIAGRAM 3236 OF 15495 ***
    // Wavefunction(s) for diagram number 3236
    // (none)
    // Amplitude(s) for diagram number 3236
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[212], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[573] -= amp_sv[0];

    // *** DIAGRAM 3237 OF 15495 ***
    // Wavefunction(s) for diagram number 3237
    // (none)
    // Amplitude(s) for diagram number 3237
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[563] -= amp_sv[0];

    // *** DIAGRAM 3238 OF 15495 ***
    // Wavefunction(s) for diagram number 3238
    // (none)
    // Amplitude(s) for diagram number 3238
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[455], w_fp[213], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[597] -= amp_sv[0];

    // *** DIAGRAM 3239 OF 15495 ***
    // Wavefunction(s) for diagram number 3239
    // (none)
    // Amplitude(s) for diagram number 3239
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[440], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[587] -= amp_sv[0];

    // *** DIAGRAM 3240 OF 15495 ***
    // Wavefunction(s) for diagram number 3240
    // (none)
    // Amplitude(s) for diagram number 3240
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[446], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[525] -= amp_sv[0];

    // *** DIAGRAM 3241 OF 15495 ***
    // Wavefunction(s) for diagram number 3241
    // (none)
    // Amplitude(s) for diagram number 3241
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= amp_sv[0];

    // *** DIAGRAM 3242 OF 15495 ***
    // Wavefunction(s) for diagram number 3242
    // (none)
    // Amplitude(s) for diagram number 3242
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[513], w_fp[211], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[549] -= amp_sv[0];

    // *** DIAGRAM 3243 OF 15495 ***
    // Wavefunction(s) for diagram number 3243
    // (none)
    // Amplitude(s) for diagram number 3243
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[539] -= amp_sv[0];

    // *** DIAGRAM 3244 OF 15495 ***
    // Wavefunction(s) for diagram number 3244
    // (none)
    // Amplitude(s) for diagram number 3244
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[513], w_fp[213], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[591] -= amp_sv[0];

    // *** DIAGRAM 3245 OF 15495 ***
    // Wavefunction(s) for diagram number 3245
    // (none)
    // Amplitude(s) for diagram number 3245
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[446], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] -= amp_sv[0];

    // *** DIAGRAM 3246 OF 15495 ***
    // Wavefunction(s) for diagram number 3246
    // (none)
    // Amplitude(s) for diagram number 3246
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] -= amp_sv[0];

    // *** DIAGRAM 3247 OF 15495 ***
    // Wavefunction(s) for diagram number 3247
    // (none)
    // Amplitude(s) for diagram number 3247
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= amp_sv[0];

    // *** DIAGRAM 3248 OF 15495 ***
    // Wavefunction(s) for diagram number 3248
    // (none)
    // Amplitude(s) for diagram number 3248
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[211], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[543] -= amp_sv[0];

    // *** DIAGRAM 3249 OF 15495 ***
    // Wavefunction(s) for diagram number 3249
    // (none)
    // Amplitude(s) for diagram number 3249
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] -= amp_sv[0];

    // *** DIAGRAM 3250 OF 15495 ***
    // Wavefunction(s) for diagram number 3250
    // (none)
    // Amplitude(s) for diagram number 3250
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[515], w_fp[212], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[567] -= amp_sv[0];

    // *** DIAGRAM 3251 OF 15495 ***
    // Wavefunction(s) for diagram number 3251
    // (none)
    // Amplitude(s) for diagram number 3251
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] -= amp_sv[0];

    // *** DIAGRAM 3252 OF 15495 ***
    // Wavefunction(s) for diagram number 3252
    // (none)
    // Amplitude(s) for diagram number 3252
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[222], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3253 OF 15495 ***
    // Wavefunction(s) for diagram number 3253
    // (none)
    // Amplitude(s) for diagram number 3253
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[213], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3254 OF 15495 ***
    // Wavefunction(s) for diagram number 3254
    // (none)
    // Amplitude(s) for diagram number 3254
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[197], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 3255 OF 15495 ***
    // Wavefunction(s) for diagram number 3255
    // (none)
    // Amplitude(s) for diagram number 3255
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[254], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3256 OF 15495 ***
    // Wavefunction(s) for diagram number 3256
    // (none)
    // Amplitude(s) for diagram number 3256
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[1], w_fp[69], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3257 OF 15495 ***
    // Wavefunction(s) for diagram number 3257
    // (none)
    // Amplitude(s) for diagram number 3257
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[66], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[66], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[66], w_fp[7], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3258 OF 15495 ***
    // Wavefunction(s) for diagram number 3258
    // (none)
    // Amplitude(s) for diagram number 3258
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[456], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3259 OF 15495 ***
    // Wavefunction(s) for diagram number 3259
    // (none)
    // Amplitude(s) for diagram number 3259
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[456], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3260 OF 15495 ***
    // Wavefunction(s) for diagram number 3260
    // (none)
    // Amplitude(s) for diagram number 3260
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[476], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3261 OF 15495 ***
    // Wavefunction(s) for diagram number 3261
    // (none)
    // Amplitude(s) for diagram number 3261
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[197], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];

    // *** DIAGRAM 3262 OF 15495 ***
    // Wavefunction(s) for diagram number 3262
    // (none)
    // Amplitude(s) for diagram number 3262
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[222], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[543] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3263 OF 15495 ***
    // Wavefunction(s) for diagram number 3263
    // (none)
    // Amplitude(s) for diagram number 3263
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[476], w_fp[69], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];

    // *** DIAGRAM 3264 OF 15495 ***
    // Wavefunction(s) for diagram number 3264
    // (none)
    // Amplitude(s) for diagram number 3264
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[213], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 3265 OF 15495 ***
    // Wavefunction(s) for diagram number 3265
    // (none)
    // Amplitude(s) for diagram number 3265
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[223], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3266 OF 15495 ***
    // Wavefunction(s) for diagram number 3266
    // (none)
    // Amplitude(s) for diagram number 3266
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[212], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3267 OF 15495 ***
    // Wavefunction(s) for diagram number 3267
    // (none)
    // Amplitude(s) for diagram number 3267
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[197], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[551] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];

    // *** DIAGRAM 3268 OF 15495 ***
    // Wavefunction(s) for diagram number 3268
    // (none)
    // Amplitude(s) for diagram number 3268
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[363], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3269 OF 15495 ***
    // Wavefunction(s) for diagram number 3269
    // (none)
    // Amplitude(s) for diagram number 3269
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[1], w_fp[103], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3270 OF 15495 ***
    // Wavefunction(s) for diagram number 3270
    // (none)
    // Amplitude(s) for diagram number 3270
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[5], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[5], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[102], w_fp[5], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3271 OF 15495 ***
    // Wavefunction(s) for diagram number 3271
    // (none)
    // Amplitude(s) for diagram number 3271
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[462], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[519] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3272 OF 15495 ***
    // Wavefunction(s) for diagram number 3272
    // (none)
    // Amplitude(s) for diagram number 3272
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[462], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3273 OF 15495 ***
    // Wavefunction(s) for diagram number 3273
    // (none)
    // Amplitude(s) for diagram number 3273
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[476], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3274 OF 15495 ***
    // Wavefunction(s) for diagram number 3274
    // (none)
    // Amplitude(s) for diagram number 3274
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[197], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[591] += amp_sv[0];

    // *** DIAGRAM 3275 OF 15495 ***
    // Wavefunction(s) for diagram number 3275
    // (none)
    // Amplitude(s) for diagram number 3275
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[223], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[549] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3276 OF 15495 ***
    // Wavefunction(s) for diagram number 3276
    // (none)
    // Amplitude(s) for diagram number 3276
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[476], w_fp[103], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];

    // *** DIAGRAM 3277 OF 15495 ***
    // Wavefunction(s) for diagram number 3277
    // (none)
    // Amplitude(s) for diagram number 3277
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[212], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 3278 OF 15495 ***
    // Wavefunction(s) for diagram number 3278
    // (none)
    // Amplitude(s) for diagram number 3278
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[211], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3279 OF 15495 ***
    // Wavefunction(s) for diagram number 3279
    // (none)
    // Amplitude(s) for diagram number 3279
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[224], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3280 OF 15495 ***
    // Wavefunction(s) for diagram number 3280
    // (none)
    // Amplitude(s) for diagram number 3280
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[197], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 3281 OF 15495 ***
    // Wavefunction(s) for diagram number 3281
    // (none)
    // Amplitude(s) for diagram number 3281
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[360], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 3282 OF 15495 ***
    // Wavefunction(s) for diagram number 3282
    // (none)
    // Amplitude(s) for diagram number 3282
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[490], w_fp[1], w_fp[124], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3283 OF 15495 ***
    // Wavefunction(s) for diagram number 3283
    // (none)
    // Amplitude(s) for diagram number 3283
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[100], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[100], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[4], w_fp[100], w_fp[490], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3284 OF 15495 ***
    // Wavefunction(s) for diagram number 3284
    // (none)
    // Amplitude(s) for diagram number 3284
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[476], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3285 OF 15495 ***
    // Wavefunction(s) for diagram number 3285
    // (none)
    // Amplitude(s) for diagram number 3285
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[197], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[521] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];

    // *** DIAGRAM 3286 OF 15495 ***
    // Wavefunction(s) for diagram number 3286
    // (none)
    // Amplitude(s) for diagram number 3286
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[224], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3287 OF 15495 ***
    // Wavefunction(s) for diagram number 3287
    // (none)
    // Amplitude(s) for diagram number 3287
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3288 OF 15495 ***
    // Wavefunction(s) for diagram number 3288
    // (none)
    // Amplitude(s) for diagram number 3288
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[461], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3289 OF 15495 ***
    // Wavefunction(s) for diagram number 3289
    // (none)
    // Amplitude(s) for diagram number 3289
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[476], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];

    // *** DIAGRAM 3290 OF 15495 ***
    // Wavefunction(s) for diagram number 3290
    // (none)
    // Amplitude(s) for diagram number 3290
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[45], w_fp[211], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 3291 OF 15495 ***
    // Wavefunction(s) for diagram number 3291
    // (none)
    // Amplitude(s) for diagram number 3291
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[197], w_fp[139], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[197], w_fp[140], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[505], w_fp[197], w_fp[141], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 3292 OF 15495 ***
    // Wavefunction(s) for diagram number 3292
    // (none)
    // Amplitude(s) for diagram number 3292
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[139], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[140], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[521] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[141], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[519] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 3293 OF 15495 ***
    // Wavefunction(s) for diagram number 3293
    // (none)
    // Amplitude(s) for diagram number 3293
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[471], w_fp[197], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[472], w_fp[197], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[473], w_fp[197], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];

    // *** DIAGRAM 3294 OF 15495 ***
    // Wavefunction(s) for diagram number 3294
    // (none)
    // Amplitude(s) for diagram number 3294
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[511], w_fp[225], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 3295 OF 15495 ***
    // Wavefunction(s) for diagram number 3295
    // (none)
    // Amplitude(s) for diagram number 3295
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[225], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[665] -= amp_sv[0];

    // *** DIAGRAM 3296 OF 15495 ***
    // Wavefunction(s) for diagram number 3296
    // (none)
    // Amplitude(s) for diagram number 3296
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[226], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 3297 OF 15495 ***
    // Wavefunction(s) for diagram number 3297
    // (none)
    // Amplitude(s) for diagram number 3297
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[226], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 3298 OF 15495 ***
    // Wavefunction(s) for diagram number 3298
    // (none)
    // Amplitude(s) for diagram number 3298
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[478], w_fp[227], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 3299 OF 15495 ***
    // Wavefunction(s) for diagram number 3299
    // (none)
    // Amplitude(s) for diagram number 3299
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[511], w_fp[227], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 3300 OF 15495 ***
    // Wavefunction(s) for diagram number 3300
    FFV1_1<W_ACCESS, CD_ACCESS>( w_fp[215], w_fp[1], COUPs[1], 1.0, cIPD[0], cIPD[1], w_fp[473] );
    // Amplitude(s) for diagram number 3300
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[440], w_fp[473], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[647] -= amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    // In CUDA, copy the local jamp to the output global-memory jamp
    constexpr int ihel0 = 0; // the allJamps buffer already points to a specific helicity _within a super-buffer for nGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcolIhelNhel( jamps, icol, ihel0, nGoodHel ) += jamp_sv[icol]; // update jamps
#else
    // In C++, copy the local jamp to the output array passed as function argument
    for( int icol = 0; icol < ncolor; icol++ )
      jamps[icol] += jamp_sv[icol]; // update jamps
#endif

#ifdef MGONGPUCPP_GPUIMPL
#ifndef MGONGPU_RDC_DIAGRAMS
    // *** STORE WAVEFUNCTIONS FOR NEXT DIAGRAM GROUPS ***
    //for( int iwf = 0; iwf < nwf; iwf++ ) storeWf( wfs, w_cx, nevt, iwf );
    storeWf( wfs, w_cx, nevt, 473 );
    storeWf( wfs, w_cx, nevt, 476 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
