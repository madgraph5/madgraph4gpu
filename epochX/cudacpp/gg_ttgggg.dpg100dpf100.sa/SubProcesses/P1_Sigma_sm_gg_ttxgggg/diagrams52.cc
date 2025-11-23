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
  diagramgroup52( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
#endif
#endif

    // *** DIAGRAM 5101 OF 15495 ***
    // Wavefunction(s) for diagram number 5101
    // (none)
    // Amplitude(s) for diagram number 5101
    FFV1_0( w_fp[94], w_fp[2], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5102 OF 15495 ***
    // Wavefunction(s) for diagram number 5102
    // (none)
    // Amplitude(s) for diagram number 5102
    FFV1_0( w_fp[77], w_fp[518], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];

    // *** DIAGRAM 5103 OF 15495 ***
    // Wavefunction(s) for diagram number 5103
    // (none)
    // Amplitude(s) for diagram number 5103
    FFV1_0( w_fp[77], w_fp[2], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5104 OF 15495 ***
    // Wavefunction(s) for diagram number 5104
    // (none)
    // Amplitude(s) for diagram number 5104
    VVVV1_0( w_fp[488], w_fp[534], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[534], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[534], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 5105 OF 15495 ***
    // Wavefunction(s) for diagram number 5105
    // (none)
    // Amplitude(s) for diagram number 5105
    VVV1_0( w_fp[534], w_fp[6], w_fp[578], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 5106 OF 15495 ***
    // Wavefunction(s) for diagram number 5106
    // (none)
    // Amplitude(s) for diagram number 5106
    VVV1_0( w_fp[534], w_fp[4], w_fp[543], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 5107 OF 15495 ***
    // Wavefunction(s) for diagram number 5107
    FFV1_1( w_fp[2], w_fp[488], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[438] );
    // Amplitude(s) for diagram number 5107
    FFV1_0( w_fp[94], w_fp[438], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 5108 OF 15495 ***
    // Wavefunction(s) for diagram number 5108
    // (none)
    // Amplitude(s) for diagram number 5108
    FFV1_0( w_fp[94], w_fp[2], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5109 OF 15495 ***
    // Wavefunction(s) for diagram number 5109
    // (none)
    // Amplitude(s) for diagram number 5109
    FFV1_0( w_fp[39], w_fp[438], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 5110 OF 15495 ***
    // Wavefunction(s) for diagram number 5110
    // (none)
    // Amplitude(s) for diagram number 5110
    FFV1_0( w_fp[39], w_fp[2], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5111 OF 15495 ***
    // Wavefunction(s) for diagram number 5111
    // (none)
    // Amplitude(s) for diagram number 5111
    VVV1_0( w_fp[535], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    VVV1_0( w_fp[553], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    VVV1_0( w_fp[495], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5112 OF 15495 ***
    // Wavefunction(s) for diagram number 5112
    // (none)
    // Amplitude(s) for diagram number 5112
    FFV1_0( w_fp[77], w_fp[2], w_fp[535], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[495], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5113 OF 15495 ***
    // Wavefunction(s) for diagram number 5113
    // (none)
    // Amplitude(s) for diagram number 5113
    VVV1_0( w_fp[545], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    VVV1_0( w_fp[528], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    VVV1_0( w_fp[447], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5114 OF 15495 ***
    // Wavefunction(s) for diagram number 5114
    // (none)
    // Amplitude(s) for diagram number 5114
    FFV1_0( w_fp[39], w_fp[2], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[39], w_fp[2], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[39], w_fp[2], w_fp[447], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5115 OF 15495 ***
    // Wavefunction(s) for diagram number 5115
    // (none)
    // Amplitude(s) for diagram number 5115
    VVV1_0( w_fp[538], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVV1_0( w_fp[527], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    VVV1_0( w_fp[533], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 5116 OF 15495 ***
    // Wavefunction(s) for diagram number 5116
    // (none)
    // Amplitude(s) for diagram number 5116
    FFV1_0( w_fp[94], w_fp[2], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5117 OF 15495 ***
    // Wavefunction(s) for diagram number 5117
    FFV1_2( w_fp[157], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[454] );
    // Amplitude(s) for diagram number 5117
    FFV1_0( w_fp[454], w_fp[161], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5118 OF 15495 ***
    // Wavefunction(s) for diagram number 5118
    // (none)
    // Amplitude(s) for diagram number 5118
    FFV1_0( w_fp[454], w_fp[163], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5119 OF 15495 ***
    // Wavefunction(s) for diagram number 5119
    FFV1_1( w_fp[156], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[515] );
    // Amplitude(s) for diagram number 5119
    FFV1_0( w_fp[39], w_fp[515], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5120 OF 15495 ***
    // Wavefunction(s) for diagram number 5120
    // (none)
    // Amplitude(s) for diagram number 5120
    FFV1_0( w_fp[77], w_fp[515], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5121 OF 15495 ***
    // Wavefunction(s) for diagram number 5121
    FFV1P0_3( w_fp[157], w_fp[156], COUPs[1], 1.0, depCoup, 0., 0., w_fp[497] );
    // Amplitude(s) for diagram number 5121
    VVV1_0( w_fp[547], w_fp[497], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5122 OF 15495 ***
    // Wavefunction(s) for diagram number 5122
    // (none)
    // Amplitude(s) for diagram number 5122
    FFV1_0( w_fp[77], w_fp[156], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];

    // *** DIAGRAM 5123 OF 15495 ***
    // Wavefunction(s) for diagram number 5123
    // (none)
    // Amplitude(s) for diagram number 5123
    FFV1_0( w_fp[157], w_fp[163], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 5124 OF 15495 ***
    // Wavefunction(s) for diagram number 5124
    // (none)
    // Amplitude(s) for diagram number 5124
    VVV1_0( w_fp[488], w_fp[497], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5125 OF 15495 ***
    // Wavefunction(s) for diagram number 5125
    // (none)
    // Amplitude(s) for diagram number 5125
    FFV1_0( w_fp[39], w_fp[156], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];

    // *** DIAGRAM 5126 OF 15495 ***
    // Wavefunction(s) for diagram number 5126
    // (none)
    // Amplitude(s) for diagram number 5126
    FFV1_0( w_fp[157], w_fp[161], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];

    // *** DIAGRAM 5127 OF 15495 ***
    // Wavefunction(s) for diagram number 5127
    // (none)
    // Amplitude(s) for diagram number 5127
    FFV1_0( w_fp[39], w_fp[163], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5128 OF 15495 ***
    // Wavefunction(s) for diagram number 5128
    // (none)
    // Amplitude(s) for diagram number 5128
    FFV1_0( w_fp[77], w_fp[161], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5129 OF 15495 ***
    // Wavefunction(s) for diagram number 5129
    // (none)
    // Amplitude(s) for diagram number 5129
    FFV1_0( w_fp[157], w_fp[156], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5130 OF 15495 ***
    // Wavefunction(s) for diagram number 5130
    // (none)
    // Amplitude(s) for diagram number 5130
    FFV1_0( w_fp[454], w_fp[156], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];

    // *** DIAGRAM 5131 OF 15495 ***
    // Wavefunction(s) for diagram number 5131
    // (none)
    // Amplitude(s) for diagram number 5131
    FFV1_0( w_fp[157], w_fp[515], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];

    // *** DIAGRAM 5132 OF 15495 ***
    // Wavefunction(s) for diagram number 5132
    // (none)
    // Amplitude(s) for diagram number 5132
    FFV1_0( w_fp[157], w_fp[156], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5133 OF 15495 ***
    // Wavefunction(s) for diagram number 5133
    // (none)
    // Amplitude(s) for diagram number 5133
    FFV1_0( w_fp[454], w_fp[211], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5134 OF 15495 ***
    // Wavefunction(s) for diagram number 5134
    // (none)
    // Amplitude(s) for diagram number 5134
    FFV1_0( w_fp[454], w_fp[213], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5135 OF 15495 ***
    // Wavefunction(s) for diagram number 5135
    FFV1_1( w_fp[197], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[514] );
    // Amplitude(s) for diagram number 5135
    FFV1_0( w_fp[94], w_fp[514], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5136 OF 15495 ***
    // Wavefunction(s) for diagram number 5136
    // (none)
    // Amplitude(s) for diagram number 5136
    FFV1_0( w_fp[77], w_fp[514], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5137 OF 15495 ***
    // Wavefunction(s) for diagram number 5137
    // (none)
    // Amplitude(s) for diagram number 5137
    VVV1_0( w_fp[577], w_fp[542], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5138 OF 15495 ***
    // Wavefunction(s) for diagram number 5138
    // (none)
    // Amplitude(s) for diagram number 5138
    FFV1_0( w_fp[77], w_fp[197], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];

    // *** DIAGRAM 5139 OF 15495 ***
    // Wavefunction(s) for diagram number 5139
    // (none)
    // Amplitude(s) for diagram number 5139
    FFV1_0( w_fp[157], w_fp[213], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];

    // *** DIAGRAM 5140 OF 15495 ***
    // Wavefunction(s) for diagram number 5140
    // (none)
    // Amplitude(s) for diagram number 5140
    VVV1_0( w_fp[488], w_fp[542], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5141 OF 15495 ***
    // Wavefunction(s) for diagram number 5141
    // (none)
    // Amplitude(s) for diagram number 5141
    FFV1_0( w_fp[94], w_fp[197], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];

    // *** DIAGRAM 5142 OF 15495 ***
    // Wavefunction(s) for diagram number 5142
    // (none)
    // Amplitude(s) for diagram number 5142
    FFV1_0( w_fp[157], w_fp[211], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];

    // *** DIAGRAM 5143 OF 15495 ***
    // Wavefunction(s) for diagram number 5143
    // (none)
    // Amplitude(s) for diagram number 5143
    FFV1_0( w_fp[94], w_fp[213], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5144 OF 15495 ***
    // Wavefunction(s) for diagram number 5144
    // (none)
    // Amplitude(s) for diagram number 5144
    FFV1_0( w_fp[77], w_fp[211], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5145 OF 15495 ***
    // Wavefunction(s) for diagram number 5145
    // (none)
    // Amplitude(s) for diagram number 5145
    FFV1_0( w_fp[157], w_fp[197], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[447], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5146 OF 15495 ***
    // Wavefunction(s) for diagram number 5146
    // (none)
    // Amplitude(s) for diagram number 5146
    FFV1_0( w_fp[454], w_fp[197], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];

    // *** DIAGRAM 5147 OF 15495 ***
    // Wavefunction(s) for diagram number 5147
    // (none)
    // Amplitude(s) for diagram number 5147
    FFV1_0( w_fp[157], w_fp[514], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];

    // *** DIAGRAM 5148 OF 15495 ***
    // Wavefunction(s) for diagram number 5148
    // (none)
    // Amplitude(s) for diagram number 5148
    FFV1_0( w_fp[157], w_fp[197], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5149 OF 15495 ***
    // Wavefunction(s) for diagram number 5149
    // (none)
    // Amplitude(s) for diagram number 5149
    FFV1_0( w_fp[454], w_fp[225], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5150 OF 15495 ***
    // Wavefunction(s) for diagram number 5150
    // (none)
    // Amplitude(s) for diagram number 5150
    FFV1_0( w_fp[454], w_fp[227], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5151 OF 15495 ***
    // Wavefunction(s) for diagram number 5151
    FFV1_1( w_fp[215], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[519] );
    // Amplitude(s) for diagram number 5151
    FFV1_0( w_fp[94], w_fp[519], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5152 OF 15495 ***
    // Wavefunction(s) for diagram number 5152
    // (none)
    // Amplitude(s) for diagram number 5152
    FFV1_0( w_fp[39], w_fp[519], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5153 OF 15495 ***
    // Wavefunction(s) for diagram number 5153
    // (none)
    // Amplitude(s) for diagram number 5153
    VVV1_0( w_fp[577], w_fp[544], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5154 OF 15495 ***
    // Wavefunction(s) for diagram number 5154
    // (none)
    // Amplitude(s) for diagram number 5154
    FFV1_0( w_fp[39], w_fp[215], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];

    // *** DIAGRAM 5155 OF 15495 ***
    // Wavefunction(s) for diagram number 5155
    // (none)
    // Amplitude(s) for diagram number 5155
    FFV1_0( w_fp[157], w_fp[227], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[701] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5156 OF 15495 ***
    // Wavefunction(s) for diagram number 5156
    // (none)
    // Amplitude(s) for diagram number 5156
    VVV1_0( w_fp[547], w_fp[544], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5157 OF 15495 ***
    // Wavefunction(s) for diagram number 5157
    // (none)
    // Amplitude(s) for diagram number 5157
    FFV1_0( w_fp[94], w_fp[215], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 5158 OF 15495 ***
    // Wavefunction(s) for diagram number 5158
    // (none)
    // Amplitude(s) for diagram number 5158
    FFV1_0( w_fp[157], w_fp[225], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];

    // *** DIAGRAM 5159 OF 15495 ***
    // Wavefunction(s) for diagram number 5159
    // (none)
    // Amplitude(s) for diagram number 5159
    FFV1_0( w_fp[94], w_fp[227], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5160 OF 15495 ***
    // Wavefunction(s) for diagram number 5160
    // (none)
    // Amplitude(s) for diagram number 5160
    FFV1_0( w_fp[39], w_fp[225], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5161 OF 15495 ***
    // Wavefunction(s) for diagram number 5161
    // (none)
    // Amplitude(s) for diagram number 5161
    FFV1_0( w_fp[157], w_fp[215], w_fp[535], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[495], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5162 OF 15495 ***
    // Wavefunction(s) for diagram number 5162
    // (none)
    // Amplitude(s) for diagram number 5162
    FFV1_0( w_fp[454], w_fp[215], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 5163 OF 15495 ***
    // Wavefunction(s) for diagram number 5163
    // (none)
    // Amplitude(s) for diagram number 5163
    FFV1_0( w_fp[157], w_fp[519], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 5164 OF 15495 ***
    // Wavefunction(s) for diagram number 5164
    // (none)
    // Amplitude(s) for diagram number 5164
    FFV1_0( w_fp[157], w_fp[215], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5165 OF 15495 ***
    // Wavefunction(s) for diagram number 5165
    // (none)
    // Amplitude(s) for diagram number 5165
    FFV1_0( w_fp[454], w_fp[118], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];

    // *** DIAGRAM 5166 OF 15495 ***
    // Wavefunction(s) for diagram number 5166
    // (none)
    // Amplitude(s) for diagram number 5166
    FFV1_0( w_fp[454], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5167 OF 15495 ***
    // Wavefunction(s) for diagram number 5167
    // (none)
    // Amplitude(s) for diagram number 5167
    FFV1_0( w_fp[110], w_fp[532], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];

    // *** DIAGRAM 5168 OF 15495 ***
    // Wavefunction(s) for diagram number 5168
    // (none)
    // Amplitude(s) for diagram number 5168
    FFV1_0( w_fp[77], w_fp[532], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];

    // *** DIAGRAM 5169 OF 15495 ***
    // Wavefunction(s) for diagram number 5169
    // (none)
    // Amplitude(s) for diagram number 5169
    FFV1_0( w_fp[157], w_fp[532], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5170 OF 15495 ***
    // Wavefunction(s) for diagram number 5170
    // (none)
    // Amplitude(s) for diagram number 5170
    VVV1_0( w_fp[571], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 5171 OF 15495 ***
    // Wavefunction(s) for diagram number 5171
    // (none)
    // Amplitude(s) for diagram number 5171
    FFV1_0( w_fp[77], w_fp[2], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5172 OF 15495 ***
    // Wavefunction(s) for diagram number 5172
    // (none)
    // Amplitude(s) for diagram number 5172
    VVV1_0( w_fp[488], w_fp[534], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];

    // *** DIAGRAM 5173 OF 15495 ***
    // Wavefunction(s) for diagram number 5173
    // (none)
    // Amplitude(s) for diagram number 5173
    FFV1_0( w_fp[110], w_fp[2], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5174 OF 15495 ***
    // Wavefunction(s) for diagram number 5174
    // (none)
    // Amplitude(s) for diagram number 5174
    FFV1_0( w_fp[157], w_fp[118], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5175 OF 15495 ***
    // Wavefunction(s) for diagram number 5175
    // (none)
    // Amplitude(s) for diagram number 5175
    VVV1_0( w_fp[479], w_fp[534], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 5176 OF 15495 ***
    // Wavefunction(s) for diagram number 5176
    // (none)
    // Amplitude(s) for diagram number 5176
    FFV1_0( w_fp[77], w_fp[118], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];

    // *** DIAGRAM 5177 OF 15495 ***
    // Wavefunction(s) for diagram number 5177
    // (none)
    // Amplitude(s) for diagram number 5177
    FFV1_0( w_fp[157], w_fp[2], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[531], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 5178 OF 15495 ***
    // Wavefunction(s) for diagram number 5178
    // (none)
    // Amplitude(s) for diagram number 5178
    FFV1_0( w_fp[454], w_fp[98], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];

    // *** DIAGRAM 5179 OF 15495 ***
    // Wavefunction(s) for diagram number 5179
    // (none)
    // Amplitude(s) for diagram number 5179
    FFV1_0( w_fp[454], w_fp[2], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5180 OF 15495 ***
    // Wavefunction(s) for diagram number 5180
    // (none)
    // Amplitude(s) for diagram number 5180
    FFV1_0( w_fp[120], w_fp[532], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];

    // *** DIAGRAM 5181 OF 15495 ***
    // Wavefunction(s) for diagram number 5181
    // (none)
    // Amplitude(s) for diagram number 5181
    FFV1_0( w_fp[39], w_fp[532], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];

    // *** DIAGRAM 5182 OF 15495 ***
    // Wavefunction(s) for diagram number 5182
    // (none)
    // Amplitude(s) for diagram number 5182
    FFV1_0( w_fp[157], w_fp[532], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5183 OF 15495 ***
    // Wavefunction(s) for diagram number 5183
    // (none)
    // Amplitude(s) for diagram number 5183
    VVV1_0( w_fp[579], w_fp[534], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];

    // *** DIAGRAM 5184 OF 15495 ***
    // Wavefunction(s) for diagram number 5184
    // (none)
    // Amplitude(s) for diagram number 5184
    FFV1_0( w_fp[39], w_fp[2], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[339] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5185 OF 15495 ***
    // Wavefunction(s) for diagram number 5185
    // (none)
    // Amplitude(s) for diagram number 5185
    VVV1_0( w_fp[547], w_fp[534], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];

    // *** DIAGRAM 5186 OF 15495 ***
    // Wavefunction(s) for diagram number 5186
    // (none)
    // Amplitude(s) for diagram number 5186
    FFV1_0( w_fp[120], w_fp[2], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[557] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5187 OF 15495 ***
    // Wavefunction(s) for diagram number 5187
    // (none)
    // Amplitude(s) for diagram number 5187
    FFV1_0( w_fp[157], w_fp[98], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5188 OF 15495 ***
    // Wavefunction(s) for diagram number 5188
    // (none)
    // Amplitude(s) for diagram number 5188
    VVV1_0( w_fp[479], w_fp[534], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];

    // *** DIAGRAM 5189 OF 15495 ***
    // Wavefunction(s) for diagram number 5189
    // (none)
    // Amplitude(s) for diagram number 5189
    FFV1_0( w_fp[39], w_fp[98], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[651] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];

    // *** DIAGRAM 5190 OF 15495 ***
    // Wavefunction(s) for diagram number 5190
    // (none)
    // Amplitude(s) for diagram number 5190
    FFV1_0( w_fp[157], w_fp[2], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[339] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[557] -= amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];

    // *** DIAGRAM 5191 OF 15495 ***
    // Wavefunction(s) for diagram number 5191
    // (none)
    // Amplitude(s) for diagram number 5191
    FFV1_0( w_fp[454], w_fp[128], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[589] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 5192 OF 15495 ***
    // Wavefunction(s) for diagram number 5192
    // (none)
    // Amplitude(s) for diagram number 5192
    FFV1_0( w_fp[454], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5193 OF 15495 ***
    // Wavefunction(s) for diagram number 5193
    // (none)
    // Amplitude(s) for diagram number 5193
    FFV1_0( w_fp[94], w_fp[532], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];

    // *** DIAGRAM 5194 OF 15495 ***
    // Wavefunction(s) for diagram number 5194
    // (none)
    // Amplitude(s) for diagram number 5194
    FFV1_0( w_fp[246], w_fp[532], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];

    // *** DIAGRAM 5195 OF 15495 ***
    // Wavefunction(s) for diagram number 5195
    // (none)
    // Amplitude(s) for diagram number 5195
    FFV1_0( w_fp[157], w_fp[532], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5196 OF 15495 ***
    // Wavefunction(s) for diagram number 5196
    // (none)
    // Amplitude(s) for diagram number 5196
    VVV1_0( w_fp[577], w_fp[534], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5197 OF 15495 ***
    // Wavefunction(s) for diagram number 5197
    // (none)
    // Amplitude(s) for diagram number 5197
    FFV1_0( w_fp[246], w_fp[2], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5198 OF 15495 ***
    // Wavefunction(s) for diagram number 5198
    // (none)
    // Amplitude(s) for diagram number 5198
    FFV1_0( w_fp[157], w_fp[128], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[581] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5199 OF 15495 ***
    // Wavefunction(s) for diagram number 5199
    // (none)
    // Amplitude(s) for diagram number 5199
    VVV1_0( w_fp[581], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[293] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5200 OF 15495 ***
    // Wavefunction(s) for diagram number 5200
    // (none)
    // Amplitude(s) for diagram number 5200
    FFV1_0( w_fp[94], w_fp[2], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[581] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 438 );
    storeWf( wfs, w_cx, nevt, 454 );
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 514 );
    storeWf( wfs, w_cx, nevt, 515 );
    storeWf( wfs, w_cx, nevt, 519 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
