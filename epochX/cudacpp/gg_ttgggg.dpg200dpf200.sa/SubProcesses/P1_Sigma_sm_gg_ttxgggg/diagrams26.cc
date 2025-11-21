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
  diagramgroup26( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
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
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 161 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 439 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 460 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 467 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 517 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
#endif
#endif

    // *** DIAGRAM 5001 OF 15495 ***
    // Wavefunction(s) for diagram number 5001
    // (none)
    // Amplitude(s) for diagram number 5001
    VVV1_0( w_fp[488], w_fp[549], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5002 OF 15495 ***
    // Wavefunction(s) for diagram number 5002
    // (none)
    // Amplitude(s) for diagram number 5002
    FFV1_0( w_fp[196], w_fp[439], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];

    // *** DIAGRAM 5003 OF 15495 ***
    // Wavefunction(s) for diagram number 5003
    // (none)
    // Amplitude(s) for diagram number 5003
    FFV1_0( w_fp[198], w_fp[259], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];

    // *** DIAGRAM 5004 OF 15495 ***
    // Wavefunction(s) for diagram number 5004
    // (none)
    // Amplitude(s) for diagram number 5004
    FFV1_0( w_fp[199], w_fp[439], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5005 OF 15495 ***
    // Wavefunction(s) for diagram number 5005
    // (none)
    // Amplitude(s) for diagram number 5005
    FFV1_0( w_fp[198], w_fp[441], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5006 OF 15495 ***
    // Wavefunction(s) for diagram number 5006
    // (none)
    // Amplitude(s) for diagram number 5006
    FFV1_0( w_fp[196], w_fp[259], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5007 OF 15495 ***
    // Wavefunction(s) for diagram number 5007
    // (none)
    // Amplitude(s) for diagram number 5007
    FFV1_0( w_fp[196], w_fp[529], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];

    // *** DIAGRAM 5008 OF 15495 ***
    // Wavefunction(s) for diagram number 5008
    // (none)
    // Amplitude(s) for diagram number 5008
    FFV1_0( w_fp[453], w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];

    // *** DIAGRAM 5009 OF 15495 ***
    // Wavefunction(s) for diagram number 5009
    VVV1P0_1( w_fp[479], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[581] );
    // Amplitude(s) for diagram number 5009
    FFV1_0( w_fp[196], w_fp[259], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5010 OF 15495 ***
    // Wavefunction(s) for diagram number 5010
    // (none)
    // Amplitude(s) for diagram number 5010
    FFV1_0( w_fp[202], w_fp[529], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5011 OF 15495 ***
    // Wavefunction(s) for diagram number 5011
    // (none)
    // Amplitude(s) for diagram number 5011
    FFV1_0( w_fp[176], w_fp[529], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5012 OF 15495 ***
    // Wavefunction(s) for diagram number 5012
    FFV1_2( w_fp[174], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[580] );
    // Amplitude(s) for diagram number 5012
    FFV1_0( w_fp[580], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5013 OF 15495 ***
    // Wavefunction(s) for diagram number 5013
    // (none)
    // Amplitude(s) for diagram number 5013
    FFV1_0( w_fp[580], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5014 OF 15495 ***
    // Wavefunction(s) for diagram number 5014
    // (none)
    // Amplitude(s) for diagram number 5014
    VVV1_0( w_fp[577], w_fp[475], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5015 OF 15495 ***
    // Wavefunction(s) for diagram number 5015
    // (none)
    // Amplitude(s) for diagram number 5015
    FFV1_0( w_fp[174], w_fp[441], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];

    // *** DIAGRAM 5016 OF 15495 ***
    // Wavefunction(s) for diagram number 5016
    // (none)
    // Amplitude(s) for diagram number 5016
    FFV1_0( w_fp[176], w_fp[259], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];

    // *** DIAGRAM 5017 OF 15495 ***
    // Wavefunction(s) for diagram number 5017
    // (none)
    // Amplitude(s) for diagram number 5017
    VVV1_0( w_fp[488], w_fp[475], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5018 OF 15495 ***
    // Wavefunction(s) for diagram number 5018
    // (none)
    // Amplitude(s) for diagram number 5018
    FFV1_0( w_fp[174], w_fp[443], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];

    // *** DIAGRAM 5019 OF 15495 ***
    // Wavefunction(s) for diagram number 5019
    // (none)
    // Amplitude(s) for diagram number 5019
    FFV1_0( w_fp[202], w_fp[259], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];

    // *** DIAGRAM 5020 OF 15495 ***
    // Wavefunction(s) for diagram number 5020
    // (none)
    // Amplitude(s) for diagram number 5020
    FFV1_0( w_fp[176], w_fp[443], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5021 OF 15495 ***
    // Wavefunction(s) for diagram number 5021
    // (none)
    // Amplitude(s) for diagram number 5021
    FFV1_0( w_fp[202], w_fp[441], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5022 OF 15495 ***
    // Wavefunction(s) for diagram number 5022
    // (none)
    // Amplitude(s) for diagram number 5022
    FFV1_0( w_fp[174], w_fp[259], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[447], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5023 OF 15495 ***
    // Wavefunction(s) for diagram number 5023
    // (none)
    // Amplitude(s) for diagram number 5023
    FFV1_0( w_fp[174], w_fp[529], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];

    // *** DIAGRAM 5024 OF 15495 ***
    // Wavefunction(s) for diagram number 5024
    // (none)
    // Amplitude(s) for diagram number 5024
    FFV1_0( w_fp[580], w_fp[259], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];

    // *** DIAGRAM 5025 OF 15495 ***
    // Wavefunction(s) for diagram number 5025
    VVV1P0_1( w_fp[479], w_fp[102], COUPs[0], 1.0, depCoup, 0., 0., w_fp[579] );
    // Amplitude(s) for diagram number 5025
    FFV1_0( w_fp[174], w_fp[259], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5026 OF 15495 ***
    // Wavefunction(s) for diagram number 5026
    // (none)
    // Amplitude(s) for diagram number 5026
    FFV1_0( w_fp[204], w_fp[529], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5027 OF 15495 ***
    // Wavefunction(s) for diagram number 5027
    // (none)
    // Amplitude(s) for diagram number 5027
    FFV1_0( w_fp[181], w_fp[529], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5028 OF 15495 ***
    // Wavefunction(s) for diagram number 5028
    FFV1_2( w_fp[179], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[536] );
    // Amplitude(s) for diagram number 5028
    FFV1_0( w_fp[536], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5029 OF 15495 ***
    // Wavefunction(s) for diagram number 5029
    // (none)
    // Amplitude(s) for diagram number 5029
    FFV1_0( w_fp[536], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5030 OF 15495 ***
    // Wavefunction(s) for diagram number 5030
    // (none)
    // Amplitude(s) for diagram number 5030
    VVV1_0( w_fp[577], w_fp[517], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5031 OF 15495 ***
    // Wavefunction(s) for diagram number 5031
    // (none)
    // Amplitude(s) for diagram number 5031
    FFV1_0( w_fp[179], w_fp[439], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[204] -= amp_sv[0];

    // *** DIAGRAM 5032 OF 15495 ***
    // Wavefunction(s) for diagram number 5032
    // (none)
    // Amplitude(s) for diagram number 5032
    FFV1_0( w_fp[181], w_fp[259], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];

    // *** DIAGRAM 5033 OF 15495 ***
    // Wavefunction(s) for diagram number 5033
    // (none)
    // Amplitude(s) for diagram number 5033
    VVV1_0( w_fp[547], w_fp[517], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5034 OF 15495 ***
    // Wavefunction(s) for diagram number 5034
    // (none)
    // Amplitude(s) for diagram number 5034
    FFV1_0( w_fp[179], w_fp[443], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];

    // *** DIAGRAM 5035 OF 15495 ***
    // Wavefunction(s) for diagram number 5035
    // (none)
    // Amplitude(s) for diagram number 5035
    FFV1_0( w_fp[204], w_fp[259], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];

    // *** DIAGRAM 5036 OF 15495 ***
    // Wavefunction(s) for diagram number 5036
    // (none)
    // Amplitude(s) for diagram number 5036
    FFV1_0( w_fp[181], w_fp[443], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5037 OF 15495 ***
    // Wavefunction(s) for diagram number 5037
    // (none)
    // Amplitude(s) for diagram number 5037
    FFV1_0( w_fp[204], w_fp[439], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5038 OF 15495 ***
    // Wavefunction(s) for diagram number 5038
    // (none)
    // Amplitude(s) for diagram number 5038
    FFV1_0( w_fp[179], w_fp[259], w_fp[535], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[495], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5039 OF 15495 ***
    // Wavefunction(s) for diagram number 5039
    // (none)
    // Amplitude(s) for diagram number 5039
    FFV1_0( w_fp[179], w_fp[529], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];

    // *** DIAGRAM 5040 OF 15495 ***
    // Wavefunction(s) for diagram number 5040
    // (none)
    // Amplitude(s) for diagram number 5040
    FFV1_0( w_fp[536], w_fp[259], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];

    // *** DIAGRAM 5041 OF 15495 ***
    // Wavefunction(s) for diagram number 5041
    VVV1P0_1( w_fp[479], w_fp[86], COUPs[0], 1.0, depCoup, 0., 0., w_fp[571] );
    // Amplitude(s) for diagram number 5041
    FFV1_0( w_fp[179], w_fp[259], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5042 OF 15495 ***
    // Wavefunction(s) for diagram number 5042
    // (none)
    // Amplitude(s) for diagram number 5042
    FFV1_0( w_fp[206], w_fp[529], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];

    // *** DIAGRAM 5043 OF 15495 ***
    // Wavefunction(s) for diagram number 5043
    // (none)
    // Amplitude(s) for diagram number 5043
    FFV1_0( w_fp[3], w_fp[529], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5044 OF 15495 ***
    // Wavefunction(s) for diagram number 5044
    // (none)
    // Amplitude(s) for diagram number 5044
    FFV1_0( w_fp[530], w_fp[460], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];

    // *** DIAGRAM 5045 OF 15495 ***
    // Wavefunction(s) for diagram number 5045
    // (none)
    // Amplitude(s) for diagram number 5045
    FFV1_0( w_fp[530], w_fp[441], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 5046 OF 15495 ***
    // Wavefunction(s) for diagram number 5046
    // (none)
    // Amplitude(s) for diagram number 5046
    FFV1_0( w_fp[530], w_fp[259], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5047 OF 15495 ***
    // Wavefunction(s) for diagram number 5047
    // (none)
    // Amplitude(s) for diagram number 5047
    VVV1_0( w_fp[571], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 5048 OF 15495 ***
    // Wavefunction(s) for diagram number 5048
    // (none)
    // Amplitude(s) for diagram number 5048
    FFV1_0( w_fp[3], w_fp[441], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5049 OF 15495 ***
    // Wavefunction(s) for diagram number 5049
    // (none)
    // Amplitude(s) for diagram number 5049
    VVV1_0( w_fp[488], w_fp[456], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];

    // *** DIAGRAM 5050 OF 15495 ***
    // Wavefunction(s) for diagram number 5050
    // (none)
    // Amplitude(s) for diagram number 5050
    FFV1_0( w_fp[3], w_fp[460], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5051 OF 15495 ***
    // Wavefunction(s) for diagram number 5051
    // (none)
    // Amplitude(s) for diagram number 5051
    FFV1_0( w_fp[206], w_fp[259], w_fp[488], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5052 OF 15495 ***
    // Wavefunction(s) for diagram number 5052
    // (none)
    // Amplitude(s) for diagram number 5052
    VVV1_0( w_fp[479], w_fp[456], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 5053 OF 15495 ***
    // Wavefunction(s) for diagram number 5053
    // (none)
    // Amplitude(s) for diagram number 5053
    FFV1_0( w_fp[206], w_fp[441], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];

    // *** DIAGRAM 5054 OF 15495 ***
    // Wavefunction(s) for diagram number 5054
    VVVV1P0_1( w_fp[479], w_fp[86], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[450] );
    VVVV3P0_1( w_fp[479], w_fp[86], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[505] );
    VVVV4P0_1( w_fp[479], w_fp[86], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[531] );
    // Amplitude(s) for diagram number 5054
    FFV1_0( w_fp[3], w_fp[259], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[259], w_fp[505], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[531], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 5055 OF 15495 ***
    // Wavefunction(s) for diagram number 5055
    // (none)
    // Amplitude(s) for diagram number 5055
    FFV1_0( w_fp[208], w_fp[529], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];

    // *** DIAGRAM 5056 OF 15495 ***
    // Wavefunction(s) for diagram number 5056
    // (none)
    // Amplitude(s) for diagram number 5056
    FFV1_0( w_fp[3], w_fp[529], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5057 OF 15495 ***
    // Wavefunction(s) for diagram number 5057
    // (none)
    // Amplitude(s) for diagram number 5057
    FFV1_0( w_fp[530], w_fp[463], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];

    // *** DIAGRAM 5058 OF 15495 ***
    // Wavefunction(s) for diagram number 5058
    // (none)
    // Amplitude(s) for diagram number 5058
    FFV1_0( w_fp[530], w_fp[439], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];

    // *** DIAGRAM 5059 OF 15495 ***
    // Wavefunction(s) for diagram number 5059
    // (none)
    // Amplitude(s) for diagram number 5059
    FFV1_0( w_fp[530], w_fp[259], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5060 OF 15495 ***
    // Wavefunction(s) for diagram number 5060
    // (none)
    // Amplitude(s) for diagram number 5060
    VVV1_0( w_fp[579], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];

    // *** DIAGRAM 5061 OF 15495 ***
    // Wavefunction(s) for diagram number 5061
    // (none)
    // Amplitude(s) for diagram number 5061
    FFV1_0( w_fp[3], w_fp[439], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5062 OF 15495 ***
    // Wavefunction(s) for diagram number 5062
    // (none)
    // Amplitude(s) for diagram number 5062
    VVV1_0( w_fp[547], w_fp[456], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];

    // *** DIAGRAM 5063 OF 15495 ***
    // Wavefunction(s) for diagram number 5063
    // (none)
    // Amplitude(s) for diagram number 5063
    FFV1_0( w_fp[3], w_fp[463], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5064 OF 15495 ***
    // Wavefunction(s) for diagram number 5064
    // (none)
    // Amplitude(s) for diagram number 5064
    FFV1_0( w_fp[208], w_fp[259], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5065 OF 15495 ***
    // Wavefunction(s) for diagram number 5065
    // (none)
    // Amplitude(s) for diagram number 5065
    VVV1_0( w_fp[479], w_fp[456], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];

    // *** DIAGRAM 5066 OF 15495 ***
    // Wavefunction(s) for diagram number 5066
    // (none)
    // Amplitude(s) for diagram number 5066
    FFV1_0( w_fp[208], w_fp[439], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];

    // *** DIAGRAM 5067 OF 15495 ***
    // Wavefunction(s) for diagram number 5067
    VVVV1P0_1( w_fp[479], w_fp[102], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[561] );
    VVVV3P0_1( w_fp[479], w_fp[102], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[584] );
    VVVV4P0_1( w_fp[479], w_fp[102], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[583] );
    // Amplitude(s) for diagram number 5067
    FFV1_0( w_fp[3], w_fp[259], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[259], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];

    // *** DIAGRAM 5068 OF 15495 ***
    // Wavefunction(s) for diagram number 5068
    // (none)
    // Amplitude(s) for diagram number 5068
    FFV1_0( w_fp[188], w_fp[529], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];

    // *** DIAGRAM 5069 OF 15495 ***
    // Wavefunction(s) for diagram number 5069
    // (none)
    // Amplitude(s) for diagram number 5069
    FFV1_0( w_fp[3], w_fp[529], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5070 OF 15495 ***
    // Wavefunction(s) for diagram number 5070
    // (none)
    // Amplitude(s) for diagram number 5070
    FFV1_0( w_fp[530], w_fp[443], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];

    // *** DIAGRAM 5071 OF 15495 ***
    // Wavefunction(s) for diagram number 5071
    // (none)
    // Amplitude(s) for diagram number 5071
    FFV1_0( w_fp[530], w_fp[467], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 5072 OF 15495 ***
    // Wavefunction(s) for diagram number 5072
    // (none)
    // Amplitude(s) for diagram number 5072
    FFV1_0( w_fp[530], w_fp[259], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5073 OF 15495 ***
    // Wavefunction(s) for diagram number 5073
    // (none)
    // Amplitude(s) for diagram number 5073
    VVV1_0( w_fp[577], w_fp[456], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5074 OF 15495 ***
    // Wavefunction(s) for diagram number 5074
    // (none)
    // Amplitude(s) for diagram number 5074
    FFV1_0( w_fp[3], w_fp[467], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5075 OF 15495 ***
    // Wavefunction(s) for diagram number 5075
    // (none)
    // Amplitude(s) for diagram number 5075
    FFV1_0( w_fp[188], w_fp[259], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5076 OF 15495 ***
    // Wavefunction(s) for diagram number 5076
    // (none)
    // Amplitude(s) for diagram number 5076
    VVV1_0( w_fp[581], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];

    // *** DIAGRAM 5077 OF 15495 ***
    // Wavefunction(s) for diagram number 5077
    // (none)
    // Amplitude(s) for diagram number 5077
    FFV1_0( w_fp[3], w_fp[443], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5078 OF 15495 ***
    // Wavefunction(s) for diagram number 5078
    // (none)
    // Amplitude(s) for diagram number 5078
    VVV1_0( w_fp[479], w_fp[456], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5079 OF 15495 ***
    // Wavefunction(s) for diagram number 5079
    // (none)
    // Amplitude(s) for diagram number 5079
    FFV1_0( w_fp[188], w_fp[443], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];

    // *** DIAGRAM 5080 OF 15495 ***
    // Wavefunction(s) for diagram number 5080
    VVVV1P0_1( w_fp[479], w_fp[4], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[582] );
    VVVV3P0_1( w_fp[479], w_fp[4], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[537] );
    VVVV4P0_1( w_fp[479], w_fp[4], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[575] );
    // Amplitude(s) for diagram number 5080
    FFV1_0( w_fp[3], w_fp[259], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 5081 OF 15495 ***
    // Wavefunction(s) for diagram number 5081
    // (none)
    // Amplitude(s) for diagram number 5081
    FFV1_0( w_fp[3], w_fp[529], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[529], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[529], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5082 OF 15495 ***
    // Wavefunction(s) for diagram number 5082
    // (none)
    // Amplitude(s) for diagram number 5082
    FFV1_0( w_fp[530], w_fp[259], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[530], w_fp[259], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[530], w_fp[259], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5083 OF 15495 ***
    // Wavefunction(s) for diagram number 5083
    VVV1P0_1( w_fp[479], w_fp[145], COUPs[0], 1.0, depCoup, 0., 0., w_fp[529] );
    VVV1P0_1( w_fp[479], w_fp[146], COUPs[0], 1.0, depCoup, 0., 0., w_fp[486] );
    VVV1P0_1( w_fp[479], w_fp[147], COUPs[0], 1.0, depCoup, 0., 0., w_fp[442] );
    // Amplitude(s) for diagram number 5083
    FFV1_0( w_fp[3], w_fp[259], w_fp[529], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[442], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 5084 OF 15495 ***
    // Wavefunction(s) for diagram number 5084
    FFV1_1( w_fp[2], w_fp[479], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[532] );
    FFV1_1( w_fp[532], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[550] );
    // Amplitude(s) for diagram number 5084
    FFV1_0( w_fp[94], w_fp[550], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5085 OF 15495 ***
    // Wavefunction(s) for diagram number 5085
    FFV1_1( w_fp[532], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[590] );
    // Amplitude(s) for diagram number 5085
    FFV1_0( w_fp[94], w_fp[590], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5086 OF 15495 ***
    // Wavefunction(s) for diagram number 5086
    FFV1_1( w_fp[532], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[589] );
    // Amplitude(s) for diagram number 5086
    FFV1_0( w_fp[39], w_fp[589], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5087 OF 15495 ***
    // Wavefunction(s) for diagram number 5087
    // (none)
    // Amplitude(s) for diagram number 5087
    FFV1_0( w_fp[39], w_fp[590], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5088 OF 15495 ***
    // Wavefunction(s) for diagram number 5088
    // (none)
    // Amplitude(s) for diagram number 5088
    FFV1_0( w_fp[77], w_fp[589], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5089 OF 15495 ***
    // Wavefunction(s) for diagram number 5089
    // (none)
    // Amplitude(s) for diagram number 5089
    FFV1_0( w_fp[77], w_fp[550], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5090 OF 15495 ***
    // Wavefunction(s) for diagram number 5090
    // (none)
    // Amplitude(s) for diagram number 5090
    VVVV1_0( w_fp[577], w_fp[534], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[577], w_fp[534], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVVV4_0( w_fp[577], w_fp[534], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 5091 OF 15495 ***
    // Wavefunction(s) for diagram number 5091
    // (none)
    // Amplitude(s) for diagram number 5091
    VVV1_0( w_fp[534], w_fp[7], w_fp[471], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5092 OF 15495 ***
    // Wavefunction(s) for diagram number 5092
    // (none)
    // Amplitude(s) for diagram number 5092
    VVV1_0( w_fp[534], w_fp[6], w_fp[539], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];

    // *** DIAGRAM 5093 OF 15495 ***
    // Wavefunction(s) for diagram number 5093
    FFV1_1( w_fp[2], w_fp[577], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[588] );
    // Amplitude(s) for diagram number 5093
    FFV1_0( w_fp[39], w_fp[588], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[293] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];

    // *** DIAGRAM 5094 OF 15495 ***
    // Wavefunction(s) for diagram number 5094
    // (none)
    // Amplitude(s) for diagram number 5094
    FFV1_0( w_fp[39], w_fp[2], w_fp[539], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[293] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[651] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5095 OF 15495 ***
    // Wavefunction(s) for diagram number 5095
    // (none)
    // Amplitude(s) for diagram number 5095
    FFV1_0( w_fp[77], w_fp[588], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];

    // *** DIAGRAM 5096 OF 15495 ***
    // Wavefunction(s) for diagram number 5096
    // (none)
    // Amplitude(s) for diagram number 5096
    FFV1_0( w_fp[77], w_fp[2], w_fp[471], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[555] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5097 OF 15495 ***
    // Wavefunction(s) for diagram number 5097
    // (none)
    // Amplitude(s) for diagram number 5097
    VVVV1_0( w_fp[547], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[547], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    VVVV4_0( w_fp[547], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5098 OF 15495 ***
    // Wavefunction(s) for diagram number 5098
    // (none)
    // Amplitude(s) for diagram number 5098
    VVV1_0( w_fp[534], w_fp[7], w_fp[435], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[63] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5099 OF 15495 ***
    // Wavefunction(s) for diagram number 5099
    // (none)
    // Amplitude(s) for diagram number 5099
    VVV1_0( w_fp[534], w_fp[4], w_fp[516], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[339] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];

    // *** DIAGRAM 5100 OF 15495 ***
    // Wavefunction(s) for diagram number 5100
    FFV1_1( w_fp[2], w_fp[547], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[518] );
    // Amplitude(s) for diagram number 5100
    FFV1_0( w_fp[94], w_fp[518], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[65] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[557] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 450 );
    storeWf( wfs, w_cx, nevt, 454 );
    storeWf( wfs, w_cx, nevt, 486 );
    storeWf( wfs, w_cx, nevt, 497 );
    storeWf( wfs, w_cx, nevt, 505 );
    storeWf( wfs, w_cx, nevt, 514 );
    storeWf( wfs, w_cx, nevt, 515 );
    storeWf( wfs, w_cx, nevt, 518 );
    storeWf( wfs, w_cx, nevt, 519 );
    storeWf( wfs, w_cx, nevt, 529 );
    storeWf( wfs, w_cx, nevt, 531 );
    storeWf( wfs, w_cx, nevt, 532 );
    storeWf( wfs, w_cx, nevt, 536 );
    storeWf( wfs, w_cx, nevt, 537 );
    storeWf( wfs, w_cx, nevt, 550 );
    storeWf( wfs, w_cx, nevt, 561 );
    storeWf( wfs, w_cx, nevt, 571 );
    storeWf( wfs, w_cx, nevt, 575 );
    storeWf( wfs, w_cx, nevt, 579 );
    storeWf( wfs, w_cx, nevt, 580 );
    storeWf( wfs, w_cx, nevt, 581 );
    storeWf( wfs, w_cx, nevt, 582 );
    storeWf( wfs, w_cx, nevt, 583 );
    storeWf( wfs, w_cx, nevt, 584 );
    storeWf( wfs, w_cx, nevt, 588 );
    storeWf( wfs, w_cx, nevt, 589 );
    storeWf( wfs, w_cx, nevt, 590 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
