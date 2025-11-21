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
  diagramgroup73( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 14 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 22 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 302 );
    retrieveWf( wfs, w_cx, nevt, 304 );
    retrieveWf( wfs, w_cx, nevt, 314 );
    retrieveWf( wfs, w_cx, nevt, 324 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 344 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 422 );
    retrieveWf( wfs, w_cx, nevt, 423 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 439 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 494 );
    retrieveWf( wfs, w_cx, nevt, 496 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 565 );
    retrieveWf( wfs, w_cx, nevt, 603 );
    retrieveWf( wfs, w_cx, nevt, 628 );
    retrieveWf( wfs, w_cx, nevt, 697 );
#endif
#endif

    // *** DIAGRAM 14401 OF 15495 ***
    // Wavefunction(s) for diagram number 14401
    // (none)
    // Amplitude(s) for diagram number 14401
    FFV1_0( w_fp[174], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[110], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];

    // *** DIAGRAM 14402 OF 15495 ***
    // Wavefunction(s) for diagram number 14402
    // (none)
    // Amplitude(s) for diagram number 14402
    FFV1_0( w_fp[221], w_fp[30], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[221], w_fp[628], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[221], w_fp[141], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14403 OF 15495 ***
    // Wavefunction(s) for diagram number 14403
    // (none)
    // Amplitude(s) for diagram number 14403
    FFV1_0( w_fp[3], w_fp[30], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[141], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 14404 OF 15495 ***
    // Wavefunction(s) for diagram number 14404
    // (none)
    // Amplitude(s) for diagram number 14404
    FFV1_0( w_fp[64], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[92], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14405 OF 15495 ***
    // Wavefunction(s) for diagram number 14405
    // (none)
    // Amplitude(s) for diagram number 14405
    FFV1_0( w_fp[64], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[92], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 14406 OF 15495 ***
    // Wavefunction(s) for diagram number 14406
    // (none)
    // Amplitude(s) for diagram number 14406
    FFV1_0( w_fp[3], w_fp[148], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[148], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[292] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[148], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];

    // *** DIAGRAM 14407 OF 15495 ***
    // Wavefunction(s) for diagram number 14407
    // (none)
    // Amplitude(s) for diagram number 14407
    FFV1_0( w_fp[221], w_fp[2], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    FFV1_0( w_fp[221], w_fp[2], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    FFV1_0( w_fp[221], w_fp[2], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 14408 OF 15495 ***
    // Wavefunction(s) for diagram number 14408
    // (none)
    // Amplitude(s) for diagram number 14408
    FFV1_0( w_fp[206], w_fp[30], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[206], w_fp[628], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[206], w_fp[141], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14409 OF 15495 ***
    // Wavefunction(s) for diagram number 14409
    // (none)
    // Amplitude(s) for diagram number 14409
    FFV1_0( w_fp[3], w_fp[30], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[141], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 14410 OF 15495 ***
    // Wavefunction(s) for diagram number 14410
    // (none)
    // Amplitude(s) for diagram number 14410
    FFV1_0( w_fp[64], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[92], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14411 OF 15495 ***
    // Wavefunction(s) for diagram number 14411
    // (none)
    // Amplitude(s) for diagram number 14411
    FFV1_0( w_fp[64], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[92], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];

    // *** DIAGRAM 14412 OF 15495 ***
    // Wavefunction(s) for diagram number 14412
    // (none)
    // Amplitude(s) for diagram number 14412
    FFV1_0( w_fp[3], w_fp[118], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[118], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[118], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];

    // *** DIAGRAM 14413 OF 15495 ***
    // Wavefunction(s) for diagram number 14413
    // (none)
    // Amplitude(s) for diagram number 14413
    FFV1_0( w_fp[206], w_fp[2], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    FFV1_0( w_fp[206], w_fp[2], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    FFV1_0( w_fp[206], w_fp[2], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];

    // *** DIAGRAM 14414 OF 15495 ***
    // Wavefunction(s) for diagram number 14414
    // (none)
    // Amplitude(s) for diagram number 14414
    FFV1_0( w_fp[184], w_fp[30], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[184], w_fp[628], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[184], w_fp[141], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14415 OF 15495 ***
    // Wavefunction(s) for diagram number 14415
    // (none)
    // Amplitude(s) for diagram number 14415
    FFV1_0( w_fp[3], w_fp[30], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[141], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 14416 OF 15495 ***
    // Wavefunction(s) for diagram number 14416
    // (none)
    // Amplitude(s) for diagram number 14416
    FFV1_0( w_fp[64], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[447] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[567] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[92], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14417 OF 15495 ***
    // Wavefunction(s) for diagram number 14417
    // (none)
    // Amplitude(s) for diagram number 14417
    FFV1_0( w_fp[64], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    FFV1_0( w_fp[92], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 14418 OF 15495 ***
    // Wavefunction(s) for diagram number 14418
    // (none)
    // Amplitude(s) for diagram number 14418
    FFV1_0( w_fp[3], w_fp[244], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[244], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[244], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];

    // *** DIAGRAM 14419 OF 15495 ***
    // Wavefunction(s) for diagram number 14419
    // (none)
    // Amplitude(s) for diagram number 14419
    FFV1_0( w_fp[184], w_fp[2], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    FFV1_0( w_fp[184], w_fp[2], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[184], w_fp[2], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];

    // *** DIAGRAM 14420 OF 15495 ***
    // Wavefunction(s) for diagram number 14420
    // (none)
    // Amplitude(s) for diagram number 14420
    FFV1_0( w_fp[3], w_fp[30], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[30], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[30], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[141], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[141], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[141], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 14421 OF 15495 ***
    // Wavefunction(s) for diagram number 14421
    // (none)
    // Amplitude(s) for diagram number 14421
    FFV1_0( w_fp[64], w_fp[2], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    FFV1_0( w_fp[92], w_fp[2], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    FFV1_0( w_fp[92], w_fp[2], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    FFV1_0( w_fp[92], w_fp[2], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];

    // *** DIAGRAM 14422 OF 15495 ***
    // Wavefunction(s) for diagram number 14422
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[92] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[135] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[134] );
    FFV1_2( w_fp[3], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[133] );
    FFV1_2( w_fp[3], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[55] );
    FFV1_2( w_fp[3], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[64] );
    // Amplitude(s) for diagram number 14422
    FFV1_0( w_fp[133], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];

    // *** DIAGRAM 14423 OF 15495 ***
    // Wavefunction(s) for diagram number 14423
    // (none)
    // Amplitude(s) for diagram number 14423
    FFV1_0( w_fp[133], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 14424 OF 15495 ***
    // Wavefunction(s) for diagram number 14424
    VVV1P0_1( w_fp[92], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[141] );
    VVV1P0_1( w_fp[135], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[628] );
    VVV1P0_1( w_fp[134], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[30] );
    // Amplitude(s) for diagram number 14424
    VVV1_0( w_fp[141], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVV1_0( w_fp[628], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    VVV1_0( w_fp[30], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14425 OF 15495 ***
    // Wavefunction(s) for diagram number 14425
    // (none)
    // Amplitude(s) for diagram number 14425
    FFV1_0( w_fp[3], w_fp[441], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14426 OF 15495 ***
    // Wavefunction(s) for diagram number 14426
    VVV1P0_1( w_fp[92], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[603] );
    VVV1P0_1( w_fp[135], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[697] );
    VVV1P0_1( w_fp[134], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[139] );
    // Amplitude(s) for diagram number 14426
    VVV1_0( w_fp[603], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    VVV1_0( w_fp[697], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    VVV1_0( w_fp[139], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14427 OF 15495 ***
    // Wavefunction(s) for diagram number 14427
    // (none)
    // Amplitude(s) for diagram number 14427
    FFV1_0( w_fp[3], w_fp[439], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[439], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[439], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14428 OF 15495 ***
    // Wavefunction(s) for diagram number 14428
    VVVV1P0_1( w_fp[92], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[115] );
    VVVV3P0_1( w_fp[92], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[490] );
    VVVV4P0_1( w_fp[92], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[20] );
    VVVV1P0_1( w_fp[135], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[147] );
    VVVV3P0_1( w_fp[135], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[87] );
    VVVV4P0_1( w_fp[135], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[145] );
    VVVV1P0_1( w_fp[134], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[21] );
    VVVV3P0_1( w_fp[134], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[22] );
    VVVV4P0_1( w_fp[134], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[68] );
    // Amplitude(s) for diagram number 14428
    FFV1_0( w_fp[3], w_fp[259], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[259], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[259], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14429 OF 15495 ***
    // Wavefunction(s) for diagram number 14429
    FFV1_1( w_fp[259], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[116] );
    FFV1_1( w_fp[259], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[110] );
    FFV1_1( w_fp[259], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[114] );
    // Amplitude(s) for diagram number 14429
    FFV1_0( w_fp[174], w_fp[116], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[110], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[114], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];

    // *** DIAGRAM 14430 OF 15495 ***
    // Wavefunction(s) for diagram number 14430
    FFV1_2( w_fp[174], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[610] );
    FFV1_2( w_fp[174], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[675] );
    FFV1_2( w_fp[174], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[627] );
    // Amplitude(s) for diagram number 14430
    FFV1_0( w_fp[610], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    FFV1_0( w_fp[675], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    FFV1_0( w_fp[627], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 14431 OF 15495 ***
    // Wavefunction(s) for diagram number 14431
    // (none)
    // Amplitude(s) for diagram number 14431
    FFV1_0( w_fp[174], w_fp[259], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14432 OF 15495 ***
    // Wavefunction(s) for diagram number 14432
    // (none)
    // Amplitude(s) for diagram number 14432
    FFV1_0( w_fp[179], w_fp[116], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[110], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[114], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];

    // *** DIAGRAM 14433 OF 15495 ***
    // Wavefunction(s) for diagram number 14433
    FFV1_2( w_fp[179], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[69] );
    FFV1_2( w_fp[179], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[497] );
    FFV1_2( w_fp[179], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[700] );
    // Amplitude(s) for diagram number 14433
    FFV1_0( w_fp[69], w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    FFV1_0( w_fp[497], w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    FFV1_0( w_fp[700], w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];

    // *** DIAGRAM 14434 OF 15495 ***
    // Wavefunction(s) for diagram number 14434
    // (none)
    // Amplitude(s) for diagram number 14434
    FFV1_0( w_fp[179], w_fp[259], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14435 OF 15495 ***
    // Wavefunction(s) for diagram number 14435
    // (none)
    // Amplitude(s) for diagram number 14435
    FFV1_0( w_fp[3], w_fp[116], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[110], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[114], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14436 OF 15495 ***
    // Wavefunction(s) for diagram number 14436
    // (none)
    // Amplitude(s) for diagram number 14436
    FFV1_0( w_fp[133], w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[64], w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14437 OF 15495 ***
    // Wavefunction(s) for diagram number 14437
    VVV1P0_1( w_fp[92], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[114] );
    VVV1P0_1( w_fp[135], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[110] );
    VVV1P0_1( w_fp[134], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[116] );
    // Amplitude(s) for diagram number 14437
    FFV1_0( w_fp[3], w_fp[259], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[110], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14438 OF 15495 ***
    // Wavefunction(s) for diagram number 14438
    FFV1_1( w_fp[2], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[151] );
    FFV1_1( w_fp[2], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[103] );
    FFV1_1( w_fp[2], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[153] );
    // Amplitude(s) for diagram number 14438
    FFV1_0( w_fp[39], w_fp[151], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    FFV1_0( w_fp[39], w_fp[103], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    FFV1_0( w_fp[39], w_fp[153], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];

    // *** DIAGRAM 14439 OF 15495 ***
    // Wavefunction(s) for diagram number 14439
    // (none)
    // Amplitude(s) for diagram number 14439
    FFV1_0( w_fp[77], w_fp[151], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    FFV1_0( w_fp[77], w_fp[103], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    FFV1_0( w_fp[77], w_fp[153], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];

    // *** DIAGRAM 14440 OF 15495 ***
    // Wavefunction(s) for diagram number 14440
    // (none)
    // Amplitude(s) for diagram number 14440
    VVV1_0( w_fp[141], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    VVV1_0( w_fp[628], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    VVV1_0( w_fp[30], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 14441 OF 15495 ***
    // Wavefunction(s) for diagram number 14441
    // (none)
    // Amplitude(s) for diagram number 14441
    FFV1_0( w_fp[77], w_fp[2], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14442 OF 15495 ***
    // Wavefunction(s) for diagram number 14442
    // (none)
    // Amplitude(s) for diagram number 14442
    VVV1_0( w_fp[603], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    VVV1_0( w_fp[697], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    VVV1_0( w_fp[139], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];

    // *** DIAGRAM 14443 OF 15495 ***
    // Wavefunction(s) for diagram number 14443
    // (none)
    // Amplitude(s) for diagram number 14443
    FFV1_0( w_fp[39], w_fp[2], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[39], w_fp[2], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[39], w_fp[2], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14444 OF 15495 ***
    // Wavefunction(s) for diagram number 14444
    // (none)
    // Amplitude(s) for diagram number 14444
    FFV1_0( w_fp[157], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[490], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 14445 OF 15495 ***
    // Wavefunction(s) for diagram number 14445
    FFV1_2( w_fp[157], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[442] );
    FFV1_2( w_fp[157], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[136] );
    FFV1_2( w_fp[157], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[124] );
    // Amplitude(s) for diagram number 14445
    FFV1_0( w_fp[442], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    FFV1_0( w_fp[136], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    FFV1_0( w_fp[124], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 14446 OF 15495 ***
    // Wavefunction(s) for diagram number 14446
    FFV1_1( w_fp[197], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[65] );
    FFV1_1( w_fp[197], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[132] );
    FFV1_1( w_fp[197], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[71] );
    // Amplitude(s) for diagram number 14446
    FFV1_0( w_fp[157], w_fp[65], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[132], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[71], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];

    // *** DIAGRAM 14447 OF 15495 ***
    // Wavefunction(s) for diagram number 14447
    // (none)
    // Amplitude(s) for diagram number 14447
    FFV1_0( w_fp[157], w_fp[197], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14448 OF 15495 ***
    // Wavefunction(s) for diagram number 14448
    // (none)
    // Amplitude(s) for diagram number 14448
    FFV1_0( w_fp[442], w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    FFV1_0( w_fp[136], w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    FFV1_0( w_fp[124], w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 14449 OF 15495 ***
    // Wavefunction(s) for diagram number 14449
    FFV1_1( w_fp[215], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[664] );
    FFV1_1( w_fp[215], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[120] );
    FFV1_1( w_fp[215], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[620] );
    // Amplitude(s) for diagram number 14449
    FFV1_0( w_fp[157], w_fp[664], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[120], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[620], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];

    // *** DIAGRAM 14450 OF 15495 ***
    // Wavefunction(s) for diagram number 14450
    // (none)
    // Amplitude(s) for diagram number 14450
    FFV1_0( w_fp[157], w_fp[215], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14451 OF 15495 ***
    // Wavefunction(s) for diagram number 14451
    // (none)
    // Amplitude(s) for diagram number 14451
    FFV1_0( w_fp[442], w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[136], w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[124], w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14452 OF 15495 ***
    // Wavefunction(s) for diagram number 14452
    // (none)
    // Amplitude(s) for diagram number 14452
    FFV1_0( w_fp[157], w_fp[151], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[103], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[153], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14453 OF 15495 ***
    // Wavefunction(s) for diagram number 14453
    // (none)
    // Amplitude(s) for diagram number 14453
    FFV1_0( w_fp[157], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[110], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 14454 OF 15495 ***
    // Wavefunction(s) for diagram number 14454
    // (none)
    // Amplitude(s) for diagram number 14454
    FFV1_0( w_fp[324], w_fp[151], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[324], w_fp[103], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[324], w_fp[153], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14455 OF 15495 ***
    // Wavefunction(s) for diagram number 14455
    // (none)
    // Amplitude(s) for diagram number 14455
    FFV1_0( w_fp[3], w_fp[151], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[103], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[153], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];

    // *** DIAGRAM 14456 OF 15495 ***
    // Wavefunction(s) for diagram number 14456
    // (none)
    // Amplitude(s) for diagram number 14456
    FFV1_0( w_fp[133], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[55], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[64], w_fp[496], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14457 OF 15495 ***
    // Wavefunction(s) for diagram number 14457
    // (none)
    // Amplitude(s) for diagram number 14457
    FFV1_0( w_fp[133], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[55], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 14458 OF 15495 ***
    // Wavefunction(s) for diagram number 14458
    // (none)
    // Amplitude(s) for diagram number 14458
    FFV1_0( w_fp[3], w_fp[496], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[496], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[496], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];

    // *** DIAGRAM 14459 OF 15495 ***
    // Wavefunction(s) for diagram number 14459
    // (none)
    // Amplitude(s) for diagram number 14459
    FFV1_0( w_fp[324], w_fp[2], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[324], w_fp[2], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    FFV1_0( w_fp[324], w_fp[2], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[609] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];

    // *** DIAGRAM 14460 OF 15495 ***
    // Wavefunction(s) for diagram number 14460
    // (none)
    // Amplitude(s) for diagram number 14460
    VVVV1_0( w_fp[92], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[92], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[92], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[135], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[135], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[135], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[134], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[134], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[134], w_fp[304], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14461 OF 15495 ***
    // Wavefunction(s) for diagram number 14461
    VVV1P0_1( w_fp[92], w_fp[304], COUPs[0], 1.0, depCoup, 0., 0., w_fp[314] );
    VVV1P0_1( w_fp[135], w_fp[304], COUPs[0], 1.0, depCoup, 0., 0., w_fp[124] );
    VVV1P0_1( w_fp[134], w_fp[304], COUPs[0], 1.0, depCoup, 0., 0., w_fp[136] );
    // Amplitude(s) for diagram number 14461
    VVV1_0( w_fp[9], w_fp[7], w_fp[314], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[136], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14462 OF 15495 ***
    // Wavefunction(s) for diagram number 14462
    VVV1P0_1( w_fp[92], w_fp[9], COUPs[0], 1.0, depCoup, 0., 0., w_fp[442] );
    VVV1P0_1( w_fp[135], w_fp[9], COUPs[0], 1.0, depCoup, 0., 0., w_fp[18] );
    VVV1P0_1( w_fp[134], w_fp[9], COUPs[0], 1.0, depCoup, 0., 0., w_fp[19] );
    // Amplitude(s) for diagram number 14462
    VVV1_0( w_fp[304], w_fp[7], w_fp[442], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[304], w_fp[7], w_fp[18], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[304], w_fp[7], w_fp[19], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14463 OF 15495 ***
    // Wavefunction(s) for diagram number 14463
    // (none)
    // Amplitude(s) for diagram number 14463
    VVV1_0( w_fp[304], w_fp[9], w_fp[603], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[304], w_fp[9], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[304], w_fp[9], w_fp[139], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14464 OF 15495 ***
    // Wavefunction(s) for diagram number 14464
    // (none)
    // Amplitude(s) for diagram number 14464
    FFV1_0( w_fp[3], w_fp[215], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 14465 OF 15495 ***
    // Wavefunction(s) for diagram number 14465
    // (none)
    // Amplitude(s) for diagram number 14465
    FFV1_0( w_fp[3], w_fp[664], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[120], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[620], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14466 OF 15495 ***
    // Wavefunction(s) for diagram number 14466
    // (none)
    // Amplitude(s) for diagram number 14466
    FFV1_0( w_fp[133], w_fp[215], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[55], w_fp[215], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[64], w_fp[215], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14467 OF 15495 ***
    // Wavefunction(s) for diagram number 14467
    // (none)
    // Amplitude(s) for diagram number 14467
    FFV1_0( w_fp[179], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[179], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];

    // *** DIAGRAM 14468 OF 15495 ***
    // Wavefunction(s) for diagram number 14468
    // (none)
    // Amplitude(s) for diagram number 14468
    FFV1_0( w_fp[179], w_fp[151], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[179], w_fp[103], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[153], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14469 OF 15495 ***
    // Wavefunction(s) for diagram number 14469
    // (none)
    // Amplitude(s) for diagram number 14469
    FFV1_0( w_fp[69], w_fp[2], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[497], w_fp[2], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[700], w_fp[2], w_fp[304], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14470 OF 15495 ***
    // Wavefunction(s) for diagram number 14470
    // (none)
    // Amplitude(s) for diagram number 14470
    FFV1_0( w_fp[344], w_fp[151], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[344], w_fp[103], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[344], w_fp[153], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14471 OF 15495 ***
    // Wavefunction(s) for diagram number 14471
    // (none)
    // Amplitude(s) for diagram number 14471
    FFV1_0( w_fp[3], w_fp[151], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[103], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[153], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];

    // *** DIAGRAM 14472 OF 15495 ***
    // Wavefunction(s) for diagram number 14472
    // (none)
    // Amplitude(s) for diagram number 14472
    FFV1_0( w_fp[133], w_fp[501], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[501], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[64], w_fp[501], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14473 OF 15495 ***
    // Wavefunction(s) for diagram number 14473
    // (none)
    // Amplitude(s) for diagram number 14473
    FFV1_0( w_fp[133], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];

    // *** DIAGRAM 14474 OF 15495 ***
    // Wavefunction(s) for diagram number 14474
    // (none)
    // Amplitude(s) for diagram number 14474
    FFV1_0( w_fp[3], w_fp[501], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[501], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[501], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];

    // *** DIAGRAM 14475 OF 15495 ***
    // Wavefunction(s) for diagram number 14475
    // (none)
    // Amplitude(s) for diagram number 14475
    FFV1_0( w_fp[344], w_fp[2], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    FFV1_0( w_fp[344], w_fp[2], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    FFV1_0( w_fp[344], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];

    // *** DIAGRAM 14476 OF 15495 ***
    // Wavefunction(s) for diagram number 14476
    // (none)
    // Amplitude(s) for diagram number 14476
    VVVV1_0( w_fp[92], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[92], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[92], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[135], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[135], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[135], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[134], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[134], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[134], w_fp[325], w_fp[9], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14477 OF 15495 ***
    // Wavefunction(s) for diagram number 14477
    VVV1P0_1( w_fp[92], w_fp[325], COUPs[0], 1.0, depCoup, 0., 0., w_fp[334] );
    VVV1P0_1( w_fp[135], w_fp[325], COUPs[0], 1.0, depCoup, 0., 0., w_fp[136] );
    VVV1P0_1( w_fp[134], w_fp[325], COUPs[0], 1.0, depCoup, 0., 0., w_fp[124] );
    // Amplitude(s) for diagram number 14477
    VVV1_0( w_fp[9], w_fp[6], w_fp[334], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[136], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14478 OF 15495 ***
    // Wavefunction(s) for diagram number 14478
    // (none)
    // Amplitude(s) for diagram number 14478
    VVV1_0( w_fp[325], w_fp[6], w_fp[442], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[6], w_fp[18], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[6], w_fp[19], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14479 OF 15495 ***
    // Wavefunction(s) for diagram number 14479
    // (none)
    // Amplitude(s) for diagram number 14479
    VVV1_0( w_fp[325], w_fp[9], w_fp[141], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[9], w_fp[628], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[325], w_fp[9], w_fp[30], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14480 OF 15495 ***
    // Wavefunction(s) for diagram number 14480
    // (none)
    // Amplitude(s) for diagram number 14480
    FFV1_0( w_fp[3], w_fp[197], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];

    // *** DIAGRAM 14481 OF 15495 ***
    // Wavefunction(s) for diagram number 14481
    // (none)
    // Amplitude(s) for diagram number 14481
    FFV1_0( w_fp[3], w_fp[65], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[132], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[71], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14482 OF 15495 ***
    // Wavefunction(s) for diagram number 14482
    // (none)
    // Amplitude(s) for diagram number 14482
    FFV1_0( w_fp[133], w_fp[197], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[197], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[64], w_fp[197], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14483 OF 15495 ***
    // Wavefunction(s) for diagram number 14483
    // (none)
    // Amplitude(s) for diagram number 14483
    FFV1_0( w_fp[174], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 14484 OF 15495 ***
    // Wavefunction(s) for diagram number 14484
    // (none)
    // Amplitude(s) for diagram number 14484
    FFV1_0( w_fp[174], w_fp[151], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[103], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[153], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14485 OF 15495 ***
    // Wavefunction(s) for diagram number 14485
    // (none)
    // Amplitude(s) for diagram number 14485
    FFV1_0( w_fp[610], w_fp[2], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[675], w_fp[2], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[627], w_fp[2], w_fp[325], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14486 OF 15495 ***
    // Wavefunction(s) for diagram number 14486
    VVV1P0_1( w_fp[92], w_fp[1], COUPs[0], 1.0, depCoup, 0., 0., w_fp[124] );
    VVV1P0_1( w_fp[135], w_fp[1], COUPs[0], 1.0, depCoup, 0., 0., w_fp[136] );
    VVV1P0_1( w_fp[134], w_fp[1], COUPs[0], 1.0, depCoup, 0., 0., w_fp[334] );
    // Amplitude(s) for diagram number 14486
    VVV1_0( w_fp[124], w_fp[14], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[136], w_fp[14], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[334], w_fp[14], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14487 OF 15495 ***
    // Wavefunction(s) for diagram number 14487
    // (none)
    // Amplitude(s) for diagram number 14487
    VVV1_0( w_fp[124], w_fp[16], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[136], w_fp[16], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[334], w_fp[16], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14488 OF 15495 ***
    // Wavefunction(s) for diagram number 14488
    // (none)
    // Amplitude(s) for diagram number 14488
    VVVV1_0( w_fp[9], w_fp[6], w_fp[7], w_fp[124], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[9], w_fp[6], w_fp[7], w_fp[124], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[9], w_fp[6], w_fp[7], w_fp[124], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[9], w_fp[6], w_fp[7], w_fp[136], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[9], w_fp[6], w_fp[7], w_fp[136], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[9], w_fp[6], w_fp[7], w_fp[136], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[9], w_fp[6], w_fp[7], w_fp[334], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[9], w_fp[6], w_fp[7], w_fp[334], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[9], w_fp[6], w_fp[7], w_fp[334], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14489 OF 15495 ***
    // Wavefunction(s) for diagram number 14489
    // (none)
    // Amplitude(s) for diagram number 14489
    VVV1_0( w_fp[141], w_fp[565], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[628], w_fp[565], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[30], w_fp[565], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14490 OF 15495 ***
    // Wavefunction(s) for diagram number 14490
    // (none)
    // Amplitude(s) for diagram number 14490
    VVV1_0( w_fp[141], w_fp[1], w_fp[16], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[628], w_fp[1], w_fp[16], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[30], w_fp[1], w_fp[16], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14491 OF 15495 ***
    // Wavefunction(s) for diagram number 14491
    // (none)
    // Amplitude(s) for diagram number 14491
    VVVV1_0( w_fp[1], w_fp[9], w_fp[7], w_fp[141], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[9], w_fp[7], w_fp[141], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[9], w_fp[7], w_fp[141], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[9], w_fp[7], w_fp[628], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[9], w_fp[7], w_fp[628], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[9], w_fp[7], w_fp[628], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[9], w_fp[7], w_fp[30], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[9], w_fp[7], w_fp[30], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[9], w_fp[7], w_fp[30], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14492 OF 15495 ***
    // Wavefunction(s) for diagram number 14492
    // (none)
    // Amplitude(s) for diagram number 14492
    VVV1_0( w_fp[603], w_fp[565], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[697], w_fp[565], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[139], w_fp[565], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14493 OF 15495 ***
    // Wavefunction(s) for diagram number 14493
    // (none)
    // Amplitude(s) for diagram number 14493
    VVV1_0( w_fp[603], w_fp[1], w_fp[14], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[697], w_fp[1], w_fp[14], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[139], w_fp[1], w_fp[14], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14494 OF 15495 ***
    // Wavefunction(s) for diagram number 14494
    // (none)
    // Amplitude(s) for diagram number 14494
    VVVV1_0( w_fp[1], w_fp[9], w_fp[6], w_fp[603], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[9], w_fp[6], w_fp[603], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[9], w_fp[6], w_fp[603], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[9], w_fp[6], w_fp[697], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[9], w_fp[6], w_fp[697], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[9], w_fp[6], w_fp[697], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[9], w_fp[6], w_fp[139], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[9], w_fp[6], w_fp[139], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[9], w_fp[6], w_fp[139], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14495 OF 15495 ***
    // Wavefunction(s) for diagram number 14495
    VVVV1P0_1( w_fp[92], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[314] );
    VVVV3P0_1( w_fp[92], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[88] );
    VVVV4P0_1( w_fp[92], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[698] );
    VVVV1P0_1( w_fp[135], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[238] );
    VVVV3P0_1( w_fp[135], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[182] );
    VVVV4P0_1( w_fp[135], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[62] );
    VVVV1P0_1( w_fp[134], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[126] );
    VVVV3P0_1( w_fp[134], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[90] );
    VVVV4P0_1( w_fp[134], w_fp[1], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[607] );
    // Amplitude(s) for diagram number 14495
    VVV1_0( w_fp[9], w_fp[7], w_fp[314], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[698], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[238], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[182], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[62], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[126], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[90], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[607], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14496 OF 15495 ***
    // Wavefunction(s) for diagram number 14496
    VVVV1P0_1( w_fp[92], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[604] );
    VVVV3P0_1( w_fp[92], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[246] );
    VVVV4P0_1( w_fp[92], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[121] );
    VVVV1P0_1( w_fp[135], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[623] );
    VVVV3P0_1( w_fp[135], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[131] );
    VVVV4P0_1( w_fp[135], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[572] );
    VVVV1P0_1( w_fp[134], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[573] );
    VVVV3P0_1( w_fp[134], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[574] );
    VVVV4P0_1( w_fp[134], w_fp[1], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[130] );
    // Amplitude(s) for diagram number 14496
    VVV1_0( w_fp[9], w_fp[6], w_fp[604], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[246], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[121], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[623], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[131], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[572], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[573], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[574], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[6], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14497 OF 15495 ***
    // Wavefunction(s) for diagram number 14497
    // (none)
    // Amplitude(s) for diagram number 14497
    VVV1_0( w_fp[1], w_fp[9], w_fp[115], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[490], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[20], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[147], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[87], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[145], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[21], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[22], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[68], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[609] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14498 OF 15495 ***
    // Wavefunction(s) for diagram number 14498
    // (none)
    // Amplitude(s) for diagram number 14498
    VVVV1_0( w_fp[92], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[92], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[92], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[135], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[135], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[135], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[134], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[134], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[134], w_fp[1], w_fp[9], w_fp[84], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14499 OF 15495 ***
    // Wavefunction(s) for diagram number 14499
    // (none)
    // Amplitude(s) for diagram number 14499
    VVV1_0( w_fp[9], w_fp[84], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[84], w_fp[136], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[84], w_fp[334], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14500 OF 15495 ***
    // Wavefunction(s) for diagram number 14500
    // (none)
    // Amplitude(s) for diagram number 14500
    VVV1_0( w_fp[1], w_fp[84], w_fp[442], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[84], w_fp[18], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[84], w_fp[19], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14501 OF 15495 ***
    // Wavefunction(s) for diagram number 14501
    // (none)
    // Amplitude(s) for diagram number 14501
    VVV1_0( w_fp[1], w_fp[9], w_fp[114], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[110], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[1], w_fp[9], w_fp[116], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14502 OF 15495 ***
    // Wavefunction(s) for diagram number 14502
    // (none)
    // Amplitude(s) for diagram number 14502
    VVV1_0( w_fp[124], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    VVV1_0( w_fp[136], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    VVV1_0( w_fp[334], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 14503 OF 15495 ***
    // Wavefunction(s) for diagram number 14503
    // (none)
    // Amplitude(s) for diagram number 14503
    FFV1_0( w_fp[3], w_fp[213], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[213], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[213], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14504 OF 15495 ***
    // Wavefunction(s) for diagram number 14504
    // (none)
    // Amplitude(s) for diagram number 14504
    FFV1_0( w_fp[133], w_fp[476], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[476], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[476], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];

    // *** DIAGRAM 14505 OF 15495 ***
    // Wavefunction(s) for diagram number 14505
    // (none)
    // Amplitude(s) for diagram number 14505
    FFV1_0( w_fp[133], w_fp[213], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[213], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[213], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];

    // *** DIAGRAM 14506 OF 15495 ***
    // Wavefunction(s) for diagram number 14506
    // (none)
    // Amplitude(s) for diagram number 14506
    FFV1_0( w_fp[3], w_fp[476], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[476], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[476], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14507 OF 15495 ***
    // Wavefunction(s) for diagram number 14507
    // (none)
    // Amplitude(s) for diagram number 14507
    VVV1_0( w_fp[603], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    VVV1_0( w_fp[697], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    VVV1_0( w_fp[139], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];

    // *** DIAGRAM 14508 OF 15495 ***
    // Wavefunction(s) for diagram number 14508
    // (none)
    // Amplitude(s) for diagram number 14508
    FFV1_0( w_fp[3], w_fp[197], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[246], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[121], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[623], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[572], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[573], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[574], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];

    // *** DIAGRAM 14509 OF 15495 ***
    // Wavefunction(s) for diagram number 14509
    // (none)
    // Amplitude(s) for diagram number 14509
    FFV1_0( w_fp[179], w_fp[197], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[197], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[197], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14510 OF 15495 ***
    // Wavefunction(s) for diagram number 14510
    // (none)
    // Amplitude(s) for diagram number 14510
    FFV1_0( w_fp[179], w_fp[65], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[132], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[71], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 14511 OF 15495 ***
    // Wavefunction(s) for diagram number 14511
    // (none)
    // Amplitude(s) for diagram number 14511
    FFV1_0( w_fp[69], w_fp[197], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    FFV1_0( w_fp[497], w_fp[197], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    FFV1_0( w_fp[700], w_fp[197], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];

    // *** DIAGRAM 14512 OF 15495 ***
    // Wavefunction(s) for diagram number 14512
    // (none)
    // Amplitude(s) for diagram number 14512
    VVV1_0( w_fp[124], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[136], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVV1_0( w_fp[334], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 14513 OF 15495 ***
    // Wavefunction(s) for diagram number 14513
    // (none)
    // Amplitude(s) for diagram number 14513
    FFV1_0( w_fp[3], w_fp[227], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[227], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14514 OF 15495 ***
    // Wavefunction(s) for diagram number 14514
    // (none)
    // Amplitude(s) for diagram number 14514
    FFV1_0( w_fp[133], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[473], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];

    // *** DIAGRAM 14515 OF 15495 ***
    // Wavefunction(s) for diagram number 14515
    // (none)
    // Amplitude(s) for diagram number 14515
    FFV1_0( w_fp[133], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[227], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

    // *** DIAGRAM 14516 OF 15495 ***
    // Wavefunction(s) for diagram number 14516
    // (none)
    // Amplitude(s) for diagram number 14516
    FFV1_0( w_fp[3], w_fp[473], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[473], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14517 OF 15495 ***
    // Wavefunction(s) for diagram number 14517
    // (none)
    // Amplitude(s) for diagram number 14517
    VVV1_0( w_fp[141], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[628], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    VVV1_0( w_fp[30], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 14518 OF 15495 ***
    // Wavefunction(s) for diagram number 14518
    // (none)
    // Amplitude(s) for diagram number 14518
    FFV1_0( w_fp[3], w_fp[215], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[698], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[215], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[607], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 14519 OF 15495 ***
    // Wavefunction(s) for diagram number 14519
    // (none)
    // Amplitude(s) for diagram number 14519
    FFV1_0( w_fp[174], w_fp[215], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[174], w_fp[215], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14520 OF 15495 ***
    // Wavefunction(s) for diagram number 14520
    // (none)
    // Amplitude(s) for diagram number 14520
    FFV1_0( w_fp[174], w_fp[664], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[120], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[620], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 14521 OF 15495 ***
    // Wavefunction(s) for diagram number 14521
    // (none)
    // Amplitude(s) for diagram number 14521
    FFV1_0( w_fp[610], w_fp[215], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    FFV1_0( w_fp[675], w_fp[215], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    FFV1_0( w_fp[627], w_fp[215], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];

    // *** DIAGRAM 14522 OF 15495 ***
    // Wavefunction(s) for diagram number 14522
    // (none)
    // Amplitude(s) for diagram number 14522
    VVV1_0( w_fp[124], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVV1_0( w_fp[136], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    VVV1_0( w_fp[334], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 14523 OF 15495 ***
    // Wavefunction(s) for diagram number 14523
    // (none)
    // Amplitude(s) for diagram number 14523
    FFV1_0( w_fp[176], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14524 OF 15495 ***
    // Wavefunction(s) for diagram number 14524
    // (none)
    // Amplitude(s) for diagram number 14524
    FFV1_0( w_fp[255], w_fp[151], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    FFV1_0( w_fp[255], w_fp[103], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    FFV1_0( w_fp[255], w_fp[153], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];

    // *** DIAGRAM 14525 OF 15495 ***
    // Wavefunction(s) for diagram number 14525
    // (none)
    // Amplitude(s) for diagram number 14525
    FFV1_0( w_fp[176], w_fp[151], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    FFV1_0( w_fp[176], w_fp[103], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    FFV1_0( w_fp[176], w_fp[153], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];

    // *** DIAGRAM 14526 OF 15495 ***
    // Wavefunction(s) for diagram number 14526
    // (none)
    // Amplitude(s) for diagram number 14526
    FFV1_0( w_fp[255], w_fp[2], w_fp[603], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[255], w_fp[2], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14527 OF 15495 ***
    // Wavefunction(s) for diagram number 14527
    // (none)
    // Amplitude(s) for diagram number 14527
    VVV1_0( w_fp[603], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0( w_fp[697], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    VVV1_0( w_fp[139], w_fp[1], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 14528 OF 15495 ***
    // Wavefunction(s) for diagram number 14528
    // (none)
    // Amplitude(s) for diagram number 14528
    FFV1_0( w_fp[174], w_fp[2], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[246], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[121], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[623], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[572], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[573], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[574], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 14529 OF 15495 ***
    // Wavefunction(s) for diagram number 14529
    // (none)
    // Amplitude(s) for diagram number 14529
    VVV1_0( w_fp[124], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVV1_0( w_fp[136], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    VVV1_0( w_fp[334], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 14530 OF 15495 ***
    // Wavefunction(s) for diagram number 14530
    // (none)
    // Amplitude(s) for diagram number 14530
    FFV1_0( w_fp[181], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14531 OF 15495 ***
    // Wavefunction(s) for diagram number 14531
    // (none)
    // Amplitude(s) for diagram number 14531
    FFV1_0( w_fp[252], w_fp[151], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    FFV1_0( w_fp[252], w_fp[103], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    FFV1_0( w_fp[252], w_fp[153], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];

    // *** DIAGRAM 14532 OF 15495 ***
    // Wavefunction(s) for diagram number 14532
    // (none)
    // Amplitude(s) for diagram number 14532
    FFV1_0( w_fp[181], w_fp[151], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    FFV1_0( w_fp[181], w_fp[103], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    FFV1_0( w_fp[181], w_fp[153], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];

    // *** DIAGRAM 14533 OF 15495 ***
    // Wavefunction(s) for diagram number 14533
    // (none)
    // Amplitude(s) for diagram number 14533
    FFV1_0( w_fp[252], w_fp[2], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[252], w_fp[2], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[252], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14534 OF 15495 ***
    // Wavefunction(s) for diagram number 14534
    // (none)
    // Amplitude(s) for diagram number 14534
    VVV1_0( w_fp[141], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    VVV1_0( w_fp[628], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    VVV1_0( w_fp[30], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 14535 OF 15495 ***
    // Wavefunction(s) for diagram number 14535
    // (none)
    // Amplitude(s) for diagram number 14535
    FFV1_0( w_fp[179], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[179], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[698], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[607], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 14536 OF 15495 ***
    // Wavefunction(s) for diagram number 14536
    // (none)
    // Amplitude(s) for diagram number 14536
    FFV1_0( w_fp[3], w_fp[128], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[128], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[128], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[579] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[718] -= amp_sv[0];

    // *** DIAGRAM 14537 OF 15495 ***
    // Wavefunction(s) for diagram number 14537
    // (none)
    // Amplitude(s) for diagram number 14537
    FFV1_0( w_fp[188], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    FFV1_0( w_fp[188], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    FFV1_0( w_fp[188], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];

    // *** DIAGRAM 14538 OF 15495 ***
    // Wavefunction(s) for diagram number 14538
    // (none)
    // Amplitude(s) for diagram number 14538
    FFV1_0( w_fp[3], w_fp[151], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[103], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[153], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];

    // *** DIAGRAM 14539 OF 15495 ***
    // Wavefunction(s) for diagram number 14539
    // (none)
    // Amplitude(s) for diagram number 14539
    FFV1_0( w_fp[188], w_fp[151], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[188], w_fp[103], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[188], w_fp[153], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14540 OF 15495 ***
    // Wavefunction(s) for diagram number 14540
    // (none)
    // Amplitude(s) for diagram number 14540
    FFV1_0( w_fp[133], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 14541 OF 15495 ***
    // Wavefunction(s) for diagram number 14541
    // (none)
    // Amplitude(s) for diagram number 14541
    FFV1_0( w_fp[133], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[64], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14542 OF 15495 ***
    // Wavefunction(s) for diagram number 14542
    // (none)
    // Amplitude(s) for diagram number 14542
    FFV1_0( w_fp[3], w_fp[151], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[151], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[151], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[103], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[103], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[103], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[153], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[153], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[153], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];

    // *** DIAGRAM 14543 OF 15495 ***
    // Wavefunction(s) for diagram number 14543
    // (none)
    // Amplitude(s) for diagram number 14543
    FFV1_0( w_fp[133], w_fp[2], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    FFV1_0( w_fp[133], w_fp[2], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    FFV1_0( w_fp[133], w_fp[2], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[585] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[525] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[705] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    FFV1_0( w_fp[64], w_fp[2], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

    // *** DIAGRAM 14544 OF 15495 ***
    // Wavefunction(s) for diagram number 14544
    // (none)
    // Amplitude(s) for diagram number 14544
    VVV1_0( w_fp[92], w_fp[111], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[92], w_fp[423], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[92], w_fp[422], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[135], w_fp[111], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[135], w_fp[423], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[585] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[135], w_fp[422], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[525] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[134], w_fp[111], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[134], w_fp[423], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[134], w_fp[422], w_fp[9], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[527] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[646] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14545 OF 15495 ***
    // Wavefunction(s) for diagram number 14545
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[134] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[135] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[92] );
    FFV1_2( w_fp[3], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[64] );
    FFV1_2( w_fp[3], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[55] );
    FFV1_2( w_fp[3], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[133] );
    // Amplitude(s) for diagram number 14545
    FFV1_0( w_fp[64], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    FFV1_0( w_fp[133], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];

    // *** DIAGRAM 14546 OF 15495 ***
    // Wavefunction(s) for diagram number 14546
    // (none)
    // Amplitude(s) for diagram number 14546
    FFV1_0( w_fp[64], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    FFV1_0( w_fp[133], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];

    // *** DIAGRAM 14547 OF 15495 ***
    // Wavefunction(s) for diagram number 14547
    VVV1P0_1( w_fp[134], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[153] );
    VVV1P0_1( w_fp[135], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[103] );
    VVV1P0_1( w_fp[92], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[151] );
    // Amplitude(s) for diagram number 14547
    VVV1_0( w_fp[153], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    VVV1_0( w_fp[103], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    VVV1_0( w_fp[151], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14548 OF 15495 ***
    // Wavefunction(s) for diagram number 14548
    // (none)
    // Amplitude(s) for diagram number 14548
    FFV1_0( w_fp[3], w_fp[441], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14549 OF 15495 ***
    // Wavefunction(s) for diagram number 14549
    VVV1P0_1( w_fp[134], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[361] );
    VVV1P0_1( w_fp[135], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[334] );
    VVV1P0_1( w_fp[92], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[136] );
    // Amplitude(s) for diagram number 14549
    VVV1_0( w_fp[361], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    VVV1_0( w_fp[334], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    VVV1_0( w_fp[136], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14550 OF 15495 ***
    // Wavefunction(s) for diagram number 14550
    // (none)
    // Amplitude(s) for diagram number 14550
    FFV1_0( w_fp[3], w_fp[436], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[436], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[436], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14551 OF 15495 ***
    // Wavefunction(s) for diagram number 14551
    VVVV1P0_1( w_fp[134], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[124] );
    VVVV3P0_1( w_fp[134], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[607] );
    VVVV4P0_1( w_fp[134], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[90] );
    VVVV1P0_1( w_fp[135], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[126] );
    VVVV3P0_1( w_fp[135], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[62] );
    VVVV4P0_1( w_fp[135], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[182] );
    VVVV1P0_1( w_fp[92], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[238] );
    VVVV3P0_1( w_fp[92], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[698] );
    VVVV4P0_1( w_fp[92], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[88] );
    // Amplitude(s) for diagram number 14551
    FFV1_0( w_fp[3], w_fp[259], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[607], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[259], w_fp[698], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[259], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14552 OF 15495 ***
    // Wavefunction(s) for diagram number 14552
    FFV1_1( w_fp[259], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[314] );
    FFV1_1( w_fp[259], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[30] );
    FFV1_1( w_fp[259], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[628] );
    // Amplitude(s) for diagram number 14552
    FFV1_0( w_fp[168], w_fp[314], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[30], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[628], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];

    // *** DIAGRAM 14553 OF 15495 ***
    // Wavefunction(s) for diagram number 14553
    FFV1_2( w_fp[168], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[141] );
    FFV1_2( w_fp[168], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[130] );
    FFV1_2( w_fp[168], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[574] );
    // Amplitude(s) for diagram number 14553
    FFV1_0( w_fp[141], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    FFV1_0( w_fp[130], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    FFV1_0( w_fp[574], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 14554 OF 15495 ***
    // Wavefunction(s) for diagram number 14554
    // (none)
    // Amplitude(s) for diagram number 14554
    FFV1_0( w_fp[168], w_fp[259], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14555 OF 15495 ***
    // Wavefunction(s) for diagram number 14555
    // (none)
    // Amplitude(s) for diagram number 14555
    FFV1_0( w_fp[179], w_fp[314], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[30], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[628], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];

    // *** DIAGRAM 14556 OF 15495 ***
    // Wavefunction(s) for diagram number 14556
    FFV1_2( w_fp[179], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[573] );
    FFV1_2( w_fp[179], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[572] );
    FFV1_2( w_fp[179], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[131] );
    // Amplitude(s) for diagram number 14556
    FFV1_0( w_fp[573], w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    FFV1_0( w_fp[572], w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    FFV1_0( w_fp[131], w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];

    // *** DIAGRAM 14557 OF 15495 ***
    // Wavefunction(s) for diagram number 14557
    // (none)
    // Amplitude(s) for diagram number 14557
    FFV1_0( w_fp[179], w_fp[259], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14558 OF 15495 ***
    // Wavefunction(s) for diagram number 14558
    // (none)
    // Amplitude(s) for diagram number 14558
    FFV1_0( w_fp[3], w_fp[314], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[30], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[628], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14559 OF 15495 ***
    // Wavefunction(s) for diagram number 14559
    // (none)
    // Amplitude(s) for diagram number 14559
    FFV1_0( w_fp[64], w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[133], w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14560 OF 15495 ***
    // Wavefunction(s) for diagram number 14560
    VVV1P0_1( w_fp[134], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[628] );
    VVV1P0_1( w_fp[135], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[30] );
    VVV1P0_1( w_fp[92], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[314] );
    // Amplitude(s) for diagram number 14560
    FFV1_0( w_fp[3], w_fp[259], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 14561 OF 15495 ***
    // Wavefunction(s) for diagram number 14561
    FFV1_1( w_fp[2], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[623] );
    FFV1_1( w_fp[2], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[121] );
    FFV1_1( w_fp[2], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[246] );
    // Amplitude(s) for diagram number 14561
    FFV1_0( w_fp[29], w_fp[623], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    FFV1_0( w_fp[29], w_fp[121], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    FFV1_0( w_fp[29], w_fp[246], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];

    // *** DIAGRAM 14562 OF 15495 ***
    // Wavefunction(s) for diagram number 14562
    // (none)
    // Amplitude(s) for diagram number 14562
    FFV1_0( w_fp[77], w_fp[623], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    FFV1_0( w_fp[77], w_fp[121], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    FFV1_0( w_fp[77], w_fp[246], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];

    // *** DIAGRAM 14563 OF 15495 ***
    // Wavefunction(s) for diagram number 14563
    // (none)
    // Amplitude(s) for diagram number 14563
    VVV1_0( w_fp[153], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVV1_0( w_fp[103], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    VVV1_0( w_fp[151], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 14564 OF 15495 ***
    // Wavefunction(s) for diagram number 14564
    // (none)
    // Amplitude(s) for diagram number 14564
    FFV1_0( w_fp[77], w_fp[2], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14565 OF 15495 ***
    // Wavefunction(s) for diagram number 14565
    // (none)
    // Amplitude(s) for diagram number 14565
    VVV1_0( w_fp[361], w_fp[534], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    VVV1_0( w_fp[334], w_fp[534], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    VVV1_0( w_fp[136], w_fp[534], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 14566 OF 15495 ***
    // Wavefunction(s) for diagram number 14566
    // (none)
    // Amplitude(s) for diagram number 14566
    FFV1_0( w_fp[29], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[29], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[29], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14567 OF 15495 ***
    // Wavefunction(s) for diagram number 14567
    // (none)
    // Amplitude(s) for diagram number 14567
    FFV1_0( w_fp[157], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[607], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[90], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[126], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[698], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 14568 OF 15495 ***
    // Wavefunction(s) for diagram number 14568
    FFV1_2( w_fp[157], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[604] );
    FFV1_2( w_fp[157], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[139] );
    FFV1_2( w_fp[157], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[697] );
    // Amplitude(s) for diagram number 14568
    FFV1_0( w_fp[604], w_fp[169], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    FFV1_0( w_fp[139], w_fp[169], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    FFV1_0( w_fp[697], w_fp[169], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 14569 OF 15495 ***
    // Wavefunction(s) for diagram number 14569
    FFV1_1( w_fp[169], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[603] );
    FFV1_1( w_fp[169], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[627] );
    FFV1_1( w_fp[169], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[675] );
    // Amplitude(s) for diagram number 14569
    FFV1_0( w_fp[157], w_fp[603], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[627], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[675], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];

    // *** DIAGRAM 14570 OF 15495 ***
    // Wavefunction(s) for diagram number 14570
    // (none)
    // Amplitude(s) for diagram number 14570
    FFV1_0( w_fp[157], w_fp[169], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14571 OF 15495 ***
    // Wavefunction(s) for diagram number 14571
    // (none)
    // Amplitude(s) for diagram number 14571
    FFV1_0( w_fp[604], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    FFV1_0( w_fp[139], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    FFV1_0( w_fp[697], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];

    // *** DIAGRAM 14572 OF 15495 ***
    // Wavefunction(s) for diagram number 14572
    FFV1_1( w_fp[215], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[610] );
    FFV1_1( w_fp[215], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[620] );
    FFV1_1( w_fp[215], w_fp[92], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[120] );
    // Amplitude(s) for diagram number 14572
    FFV1_0( w_fp[157], w_fp[610], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[620], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[120], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 14573 OF 15495 ***
    // Wavefunction(s) for diagram number 14573
    // (none)
    // Amplitude(s) for diagram number 14573
    FFV1_0( w_fp[157], w_fp[215], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14574 OF 15495 ***
    // Wavefunction(s) for diagram number 14574
    // (none)
    // Amplitude(s) for diagram number 14574
    FFV1_0( w_fp[604], w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[139], w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[697], w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14575 OF 15495 ***
    // Wavefunction(s) for diagram number 14575
    // (none)
    // Amplitude(s) for diagram number 14575
    FFV1_0( w_fp[157], w_fp[623], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[121], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[246], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14576 OF 15495 ***
    // Wavefunction(s) for diagram number 14576
    // (none)
    // Amplitude(s) for diagram number 14576
    FFV1_0( w_fp[157], w_fp[2], w_fp[628], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[314], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 14577 OF 15495 ***
    // Wavefunction(s) for diagram number 14577
    // (none)
    // Amplitude(s) for diagram number 14577
    FFV1_0( w_fp[302], w_fp[623], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[302], w_fp[121], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[302], w_fp[246], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14578 OF 15495 ***
    // Wavefunction(s) for diagram number 14578
    // (none)
    // Amplitude(s) for diagram number 14578
    FFV1_0( w_fp[3], w_fp[623], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[121], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[246], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];

    // *** DIAGRAM 14579 OF 15495 ***
    // Wavefunction(s) for diagram number 14579
    // (none)
    // Amplitude(s) for diagram number 14579
    FFV1_0( w_fp[64], w_fp[494], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[494], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[133], w_fp[494], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14580 OF 15495 ***
    // Wavefunction(s) for diagram number 14580
    // (none)
    // Amplitude(s) for diagram number 14580
    FFV1_0( w_fp[64], w_fp[2], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    FFV1_0( w_fp[133], w_fp[2], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 14581 OF 15495 ***
    // Wavefunction(s) for diagram number 14581
    // (none)
    // Amplitude(s) for diagram number 14581
    FFV1_0( w_fp[3], w_fp[494], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[494], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[494], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 14582 OF 15495 ***
    // Wavefunction(s) for diagram number 14582
    // (none)
    // Amplitude(s) for diagram number 14582
    FFV1_0( w_fp[302], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    FFV1_0( w_fp[302], w_fp[2], w_fp[334], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    FFV1_0( w_fp[302], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 14583 OF 15495 ***
    // Wavefunction(s) for diagram number 14583
    // (none)
    // Amplitude(s) for diagram number 14583
    VVVV1_0( w_fp[134], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[134], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[134], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[135], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[135], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[135], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[92], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[92], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[92], w_fp[279], w_fp[9], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14584 OF 15495 ***
    // Wavefunction(s) for diagram number 14584
    VVV1P0_1( w_fp[134], w_fp[279], COUPs[0], 1.0, depCoup, 0., 0., w_fp[292] );
    VVV1P0_1( w_fp[135], w_fp[279], COUPs[0], 1.0, depCoup, 0., 0., w_fp[697] );
    VVV1P0_1( w_fp[92], w_fp[279], COUPs[0], 1.0, depCoup, 0., 0., w_fp[139] );
    // Amplitude(s) for diagram number 14584
    VVV1_0( w_fp[9], w_fp[7], w_fp[292], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[7], w_fp[139], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14585 OF 15495 ***
    // Wavefunction(s) for diagram number 14585
    VVV1P0_1( w_fp[134], w_fp[9], COUPs[0], 1.0, depCoup, 0., 0., w_fp[604] );
    VVV1P0_1( w_fp[135], w_fp[9], COUPs[0], 1.0, depCoup, 0., 0., w_fp[664] );
    VVV1P0_1( w_fp[92], w_fp[9], COUPs[0], 1.0, depCoup, 0., 0., w_fp[700] );
    // Amplitude(s) for diagram number 14585
    VVV1_0( w_fp[279], w_fp[7], w_fp[604], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[279], w_fp[7], w_fp[664], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[279], w_fp[7], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14586 OF 15495 ***
    // Wavefunction(s) for diagram number 14586
    // (none)
    // Amplitude(s) for diagram number 14586
    VVV1_0( w_fp[279], w_fp[9], w_fp[361], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[279], w_fp[9], w_fp[334], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[279], w_fp[9], w_fp[136], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14587 OF 15495 ***
    // Wavefunction(s) for diagram number 14587
    // (none)
    // Amplitude(s) for diagram number 14587
    FFV1_0( w_fp[3], w_fp[215], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 14588 OF 15495 ***
    // Wavefunction(s) for diagram number 14588
    // (none)
    // Amplitude(s) for diagram number 14588
    FFV1_0( w_fp[3], w_fp[610], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[620], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[120], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[611] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14589 OF 15495 ***
    // Wavefunction(s) for diagram number 14589
    // (none)
    // Amplitude(s) for diagram number 14589
    FFV1_0( w_fp[64], w_fp[215], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[215], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[133], w_fp[215], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14590 OF 15495 ***
    // Wavefunction(s) for diagram number 14590
    // (none)
    // Amplitude(s) for diagram number 14590
    FFV1_0( w_fp[179], w_fp[2], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[697], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 14591 OF 15495 ***
    // Wavefunction(s) for diagram number 14591
    // (none)
    // Amplitude(s) for diagram number 14591
    FFV1_0( w_fp[179], w_fp[623], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[121], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[246], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14592 OF 15495 ***
    // Wavefunction(s) for diagram number 14592
    // (none)
    // Amplitude(s) for diagram number 14592
    FFV1_0( w_fp[573], w_fp[2], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[572], w_fp[2], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[131], w_fp[2], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14593 OF 15495 ***
    // Wavefunction(s) for diagram number 14593
    // (none)
    // Amplitude(s) for diagram number 14593
    FFV1_0( w_fp[344], w_fp[623], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[344], w_fp[121], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[344], w_fp[246], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14594 OF 15495 ***
    // Wavefunction(s) for diagram number 14594
    // (none)
    // Amplitude(s) for diagram number 14594
    FFV1_0( w_fp[3], w_fp[623], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[121], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[246], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];

    // *** DIAGRAM 14595 OF 15495 ***
    // Wavefunction(s) for diagram number 14595
    // (none)
    // Amplitude(s) for diagram number 14595
    FFV1_0( w_fp[64], w_fp[501], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[55], w_fp[501], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[133], w_fp[501], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14596 OF 15495 ***
    // Wavefunction(s) for diagram number 14596
    // (none)
    // Amplitude(s) for diagram number 14596
    FFV1_0( w_fp[64], w_fp[2], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    FFV1_0( w_fp[55], w_fp[2], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    FFV1_0( w_fp[133], w_fp[2], w_fp[333], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];

    // *** DIAGRAM 14597 OF 15495 ***
    // Wavefunction(s) for diagram number 14597
    // (none)
    // Amplitude(s) for diagram number 14597
    FFV1_0( w_fp[3], w_fp[501], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[501], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[501], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];

    // *** DIAGRAM 14598 OF 15495 ***
    // Wavefunction(s) for diagram number 14598
    // (none)
    // Amplitude(s) for diagram number 14598
    FFV1_0( w_fp[344], w_fp[2], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    FFV1_0( w_fp[344], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    FFV1_0( w_fp[344], w_fp[2], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];

    // *** DIAGRAM 14599 OF 15495 ***
    // Wavefunction(s) for diagram number 14599
    // (none)
    // Amplitude(s) for diagram number 14599
    VVVV1_0( w_fp[134], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[134], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[134], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[135], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[135], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[135], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV1_0( w_fp[92], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[92], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[92], w_fp[325], w_fp[9], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 14600 OF 15495 ***
    // Wavefunction(s) for diagram number 14600
    VVV1P0_1( w_fp[134], w_fp[325], COUPs[0], 1.0, depCoup, 0., 0., w_fp[333] );
    VVV1P0_1( w_fp[135], w_fp[325], COUPs[0], 1.0, depCoup, 0., 0., w_fp[139] );
    VVV1P0_1( w_fp[92], w_fp[325], COUPs[0], 1.0, depCoup, 0., 0., w_fp[697] );
    // Amplitude(s) for diagram number 14600
    VVV1_0( w_fp[9], w_fp[5], w_fp[333], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[5], w_fp[139], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[9], w_fp[5], w_fp[697], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[445] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 18 );
    storeWf( wfs, w_cx, nevt, 19 );
    storeWf( wfs, w_cx, nevt, 20 );
    storeWf( wfs, w_cx, nevt, 21 );
    storeWf( wfs, w_cx, nevt, 22 );
    storeWf( wfs, w_cx, nevt, 30 );
    storeWf( wfs, w_cx, nevt, 55 );
    storeWf( wfs, w_cx, nevt, 62 );
    storeWf( wfs, w_cx, nevt, 64 );
    storeWf( wfs, w_cx, nevt, 65 );
    storeWf( wfs, w_cx, nevt, 68 );
    storeWf( wfs, w_cx, nevt, 69 );
    storeWf( wfs, w_cx, nevt, 71 );
    storeWf( wfs, w_cx, nevt, 87 );
    storeWf( wfs, w_cx, nevt, 88 );
    storeWf( wfs, w_cx, nevt, 90 );
    storeWf( wfs, w_cx, nevt, 92 );
    storeWf( wfs, w_cx, nevt, 103 );
    storeWf( wfs, w_cx, nevt, 110 );
    storeWf( wfs, w_cx, nevt, 114 );
    storeWf( wfs, w_cx, nevt, 115 );
    storeWf( wfs, w_cx, nevt, 116 );
    storeWf( wfs, w_cx, nevt, 120 );
    storeWf( wfs, w_cx, nevt, 121 );
    storeWf( wfs, w_cx, nevt, 124 );
    storeWf( wfs, w_cx, nevt, 126 );
    storeWf( wfs, w_cx, nevt, 130 );
    storeWf( wfs, w_cx, nevt, 131 );
    storeWf( wfs, w_cx, nevt, 132 );
    storeWf( wfs, w_cx, nevt, 133 );
    storeWf( wfs, w_cx, nevt, 134 );
    storeWf( wfs, w_cx, nevt, 135 );
    storeWf( wfs, w_cx, nevt, 136 );
    storeWf( wfs, w_cx, nevt, 139 );
    storeWf( wfs, w_cx, nevt, 141 );
    storeWf( wfs, w_cx, nevt, 145 );
    storeWf( wfs, w_cx, nevt, 147 );
    storeWf( wfs, w_cx, nevt, 151 );
    storeWf( wfs, w_cx, nevt, 153 );
    storeWf( wfs, w_cx, nevt, 182 );
    storeWf( wfs, w_cx, nevt, 238 );
    storeWf( wfs, w_cx, nevt, 246 );
    storeWf( wfs, w_cx, nevt, 292 );
    storeWf( wfs, w_cx, nevt, 314 );
    storeWf( wfs, w_cx, nevt, 333 );
    storeWf( wfs, w_cx, nevt, 334 );
    storeWf( wfs, w_cx, nevt, 361 );
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 490 );
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 572 );
    storeWf( wfs, w_cx, nevt, 573 );
    storeWf( wfs, w_cx, nevt, 574 );
    storeWf( wfs, w_cx, nevt, 603 );
    storeWf( wfs, w_cx, nevt, 604 );
    storeWf( wfs, w_cx, nevt, 607 );
    storeWf( wfs, w_cx, nevt, 610 );
    storeWf( wfs, w_cx, nevt, 620 );
    storeWf( wfs, w_cx, nevt, 623 );
    storeWf( wfs, w_cx, nevt, 627 );
    storeWf( wfs, w_cx, nevt, 628 );
    storeWf( wfs, w_cx, nevt, 664 );
    storeWf( wfs, w_cx, nevt, 675 );
    storeWf( wfs, w_cx, nevt, 697 );
    storeWf( wfs, w_cx, nevt, 698 );
    storeWf( wfs, w_cx, nevt, 700 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
