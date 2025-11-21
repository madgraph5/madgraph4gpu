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
  diagramgroup3( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 48 );
    retrieveWf( wfs, w_cx, nevt, 49 );
    retrieveWf( wfs, w_cx, nevt, 50 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 56 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 83 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 89 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 93 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 105 );
    retrieveWf( wfs, w_cx, nevt, 107 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 112 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 142 );
    retrieveWf( wfs, w_cx, nevt, 143 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 149 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 165 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 172 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 192 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 195 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 202 );
#endif
#endif

    // *** DIAGRAM 401 OF 15495 ***
    // Wavefunction(s) for diagram number 401
    // (none)
    // Amplitude(s) for diagram number 401
    FFV1_0( w_fp[174], w_fp[190], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];

    // *** DIAGRAM 402 OF 15495 ***
    // Wavefunction(s) for diagram number 402
    // (none)
    // Amplitude(s) for diagram number 402
    FFV1_0( w_fp[202], w_fp[169], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];

    // *** DIAGRAM 403 OF 15495 ***
    // Wavefunction(s) for diagram number 403
    // (none)
    // Amplitude(s) for diagram number 403
    FFV1_0( w_fp[176], w_fp[190], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 404 OF 15495 ***
    // Wavefunction(s) for diagram number 404
    // (none)
    // Amplitude(s) for diagram number 404
    FFV1_0( w_fp[202], w_fp[193], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 405 OF 15495 ***
    // Wavefunction(s) for diagram number 405
    // (none)
    // Amplitude(s) for diagram number 405
    FFV1_0( w_fp[174], w_fp[169], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 406 OF 15495 ***
    // Wavefunction(s) for diagram number 406
    // (none)
    // Amplitude(s) for diagram number 406
    FFV1_0( w_fp[174], w_fp[197], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];

    // *** DIAGRAM 407 OF 15495 ***
    // Wavefunction(s) for diagram number 407
    // (none)
    // Amplitude(s) for diagram number 407
    FFV1_0( w_fp[177], w_fp[169], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 408 OF 15495 ***
    // Wavefunction(s) for diagram number 408
    // (none)
    // Amplitude(s) for diagram number 408
    FFV1_0( w_fp[174], w_fp[169], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 409 OF 15495 ***
    // Wavefunction(s) for diagram number 409
    FFV1_2( w_fp[179], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[204] );
    // Amplitude(s) for diagram number 409
    FFV1_0( w_fp[204], w_fp[197], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 410 OF 15495 ***
    // Wavefunction(s) for diagram number 410
    // (none)
    // Amplitude(s) for diagram number 410
    FFV1_0( w_fp[181], w_fp[197], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 411 OF 15495 ***
    // Wavefunction(s) for diagram number 411
    // (none)
    // Amplitude(s) for diagram number 411
    FFV1_0( w_fp[182], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 412 OF 15495 ***
    // Wavefunction(s) for diagram number 412
    // (none)
    // Amplitude(s) for diagram number 412
    FFV1_0( w_fp[182], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 413 OF 15495 ***
    // Wavefunction(s) for diagram number 413
    FFV1P0_3( w_fp[179], w_fp[169], COUPs[1], 1.0, depCoup, 0., 0., w_fp[205] );
    // Amplitude(s) for diagram number 413
    VVV1_0( w_fp[10], w_fp[205], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 414 OF 15495 ***
    // Wavefunction(s) for diagram number 414
    // (none)
    // Amplitude(s) for diagram number 414
    FFV1_0( w_fp[179], w_fp[191], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];

    // *** DIAGRAM 415 OF 15495 ***
    // Wavefunction(s) for diagram number 415
    // (none)
    // Amplitude(s) for diagram number 415
    FFV1_0( w_fp[181], w_fp[169], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];

    // *** DIAGRAM 416 OF 15495 ***
    // Wavefunction(s) for diagram number 416
    // (none)
    // Amplitude(s) for diagram number 416
    VVV1_0( w_fp[37], w_fp[205], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 417 OF 15495 ***
    // Wavefunction(s) for diagram number 417
    // (none)
    // Amplitude(s) for diagram number 417
    FFV1_0( w_fp[179], w_fp[190], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];

    // *** DIAGRAM 418 OF 15495 ***
    // Wavefunction(s) for diagram number 418
    // (none)
    // Amplitude(s) for diagram number 418
    FFV1_0( w_fp[204], w_fp[169], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];

    // *** DIAGRAM 419 OF 15495 ***
    // Wavefunction(s) for diagram number 419
    // (none)
    // Amplitude(s) for diagram number 419
    FFV1_0( w_fp[181], w_fp[190], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 420 OF 15495 ***
    // Wavefunction(s) for diagram number 420
    // (none)
    // Amplitude(s) for diagram number 420
    FFV1_0( w_fp[204], w_fp[191], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 421 OF 15495 ***
    // Wavefunction(s) for diagram number 421
    // (none)
    // Amplitude(s) for diagram number 421
    FFV1_0( w_fp[179], w_fp[169], w_fp[51], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[169], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[169], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 422 OF 15495 ***
    // Wavefunction(s) for diagram number 422
    // (none)
    // Amplitude(s) for diagram number 422
    FFV1_0( w_fp[179], w_fp[197], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];

    // *** DIAGRAM 423 OF 15495 ***
    // Wavefunction(s) for diagram number 423
    // (none)
    // Amplitude(s) for diagram number 423
    FFV1_0( w_fp[182], w_fp[169], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];

    // *** DIAGRAM 424 OF 15495 ***
    // Wavefunction(s) for diagram number 424
    // (none)
    // Amplitude(s) for diagram number 424
    FFV1_0( w_fp[179], w_fp[169], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 425 OF 15495 ***
    // Wavefunction(s) for diagram number 425
    FFV1_2( w_fp[3], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[206] );
    // Amplitude(s) for diagram number 425
    FFV1_0( w_fp[206], w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 426 OF 15495 ***
    // Wavefunction(s) for diagram number 426
    // (none)
    // Amplitude(s) for diagram number 426
    FFV1_0( w_fp[3], w_fp[197], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 427 OF 15495 ***
    // Wavefunction(s) for diagram number 427
    FFV1_1( w_fp[169], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[207] );
    // Amplitude(s) for diagram number 427
    FFV1_0( w_fp[157], w_fp[207], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];

    // *** DIAGRAM 428 OF 15495 ***
    // Wavefunction(s) for diagram number 428
    // (none)
    // Amplitude(s) for diagram number 428
    FFV1_0( w_fp[157], w_fp[193], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 429 OF 15495 ***
    // Wavefunction(s) for diagram number 429
    // (none)
    // Amplitude(s) for diagram number 429
    FFV1_0( w_fp[157], w_fp[169], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 430 OF 15495 ***
    // Wavefunction(s) for diagram number 430
    // (none)
    // Amplitude(s) for diagram number 430
    VVV1_0( w_fp[89], w_fp[194], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 431 OF 15495 ***
    // Wavefunction(s) for diagram number 431
    // (none)
    // Amplitude(s) for diagram number 431
    FFV1_0( w_fp[3], w_fp[193], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 432 OF 15495 ***
    // Wavefunction(s) for diagram number 432
    // (none)
    // Amplitude(s) for diagram number 432
    VVV1_0( w_fp[44], w_fp[194], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];

    // *** DIAGRAM 433 OF 15495 ***
    // Wavefunction(s) for diagram number 433
    // (none)
    // Amplitude(s) for diagram number 433
    FFV1_0( w_fp[3], w_fp[207], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 434 OF 15495 ***
    // Wavefunction(s) for diagram number 434
    // (none)
    // Amplitude(s) for diagram number 434
    FFV1_0( w_fp[206], w_fp[169], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 435 OF 15495 ***
    // Wavefunction(s) for diagram number 435
    // (none)
    // Amplitude(s) for diagram number 435
    VVV1_0( w_fp[8], w_fp[194], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 436 OF 15495 ***
    // Wavefunction(s) for diagram number 436
    // (none)
    // Amplitude(s) for diagram number 436
    FFV1_0( w_fp[206], w_fp[193], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];

    // *** DIAGRAM 437 OF 15495 ***
    // Wavefunction(s) for diagram number 437
    // (none)
    // Amplitude(s) for diagram number 437
    FFV1_0( w_fp[3], w_fp[169], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[98], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[99], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 438 OF 15495 ***
    // Wavefunction(s) for diagram number 438
    FFV1_2( w_fp[3], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[208] );
    // Amplitude(s) for diagram number 438
    FFV1_0( w_fp[208], w_fp[197], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];

    // *** DIAGRAM 439 OF 15495 ***
    // Wavefunction(s) for diagram number 439
    // (none)
    // Amplitude(s) for diagram number 439
    FFV1_0( w_fp[3], w_fp[197], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 440 OF 15495 ***
    // Wavefunction(s) for diagram number 440
    FFV1_1( w_fp[169], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[209] );
    // Amplitude(s) for diagram number 440
    FFV1_0( w_fp[157], w_fp[209], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];

    // *** DIAGRAM 441 OF 15495 ***
    // Wavefunction(s) for diagram number 441
    // (none)
    // Amplitude(s) for diagram number 441
    FFV1_0( w_fp[157], w_fp[191], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];

    // *** DIAGRAM 442 OF 15495 ***
    // Wavefunction(s) for diagram number 442
    // (none)
    // Amplitude(s) for diagram number 442
    FFV1_0( w_fp[157], w_fp[169], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 443 OF 15495 ***
    // Wavefunction(s) for diagram number 443
    // (none)
    // Amplitude(s) for diagram number 443
    VVV1_0( w_fp[105], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];

    // *** DIAGRAM 444 OF 15495 ***
    // Wavefunction(s) for diagram number 444
    // (none)
    // Amplitude(s) for diagram number 444
    FFV1_0( w_fp[3], w_fp[191], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 445 OF 15495 ***
    // Wavefunction(s) for diagram number 445
    // (none)
    // Amplitude(s) for diagram number 445
    VVV1_0( w_fp[37], w_fp[194], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];

    // *** DIAGRAM 446 OF 15495 ***
    // Wavefunction(s) for diagram number 446
    // (none)
    // Amplitude(s) for diagram number 446
    FFV1_0( w_fp[3], w_fp[209], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 447 OF 15495 ***
    // Wavefunction(s) for diagram number 447
    // (none)
    // Amplitude(s) for diagram number 447
    FFV1_0( w_fp[208], w_fp[169], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 448 OF 15495 ***
    // Wavefunction(s) for diagram number 448
    // (none)
    // Amplitude(s) for diagram number 448
    VVV1_0( w_fp[8], w_fp[194], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];

    // *** DIAGRAM 449 OF 15495 ***
    // Wavefunction(s) for diagram number 449
    // (none)
    // Amplitude(s) for diagram number 449
    FFV1_0( w_fp[208], w_fp[191], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];

    // *** DIAGRAM 450 OF 15495 ***
    // Wavefunction(s) for diagram number 450
    // (none)
    // Amplitude(s) for diagram number 450
    FFV1_0( w_fp[3], w_fp[169], w_fp[110], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[112], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];

    // *** DIAGRAM 451 OF 15495 ***
    // Wavefunction(s) for diagram number 451
    // (none)
    // Amplitude(s) for diagram number 451
    FFV1_0( w_fp[188], w_fp[197], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];

    // *** DIAGRAM 452 OF 15495 ***
    // Wavefunction(s) for diagram number 452
    // (none)
    // Amplitude(s) for diagram number 452
    FFV1_0( w_fp[3], w_fp[197], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 453 OF 15495 ***
    // Wavefunction(s) for diagram number 453
    // (none)
    // Amplitude(s) for diagram number 453
    FFV1_0( w_fp[157], w_fp[190], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];

    // *** DIAGRAM 454 OF 15495 ***
    // Wavefunction(s) for diagram number 454
    FFV1_1( w_fp[169], w_fp[84], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[210] );
    // Amplitude(s) for diagram number 454
    FFV1_0( w_fp[157], w_fp[210], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 455 OF 15495 ***
    // Wavefunction(s) for diagram number 455
    // (none)
    // Amplitude(s) for diagram number 455
    FFV1_0( w_fp[157], w_fp[169], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 456 OF 15495 ***
    // Wavefunction(s) for diagram number 456
    // (none)
    // Amplitude(s) for diagram number 456
    VVV1_0( w_fp[10], w_fp[194], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 457 OF 15495 ***
    // Wavefunction(s) for diagram number 457
    // (none)
    // Amplitude(s) for diagram number 457
    FFV1_0( w_fp[3], w_fp[210], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 458 OF 15495 ***
    // Wavefunction(s) for diagram number 458
    // (none)
    // Amplitude(s) for diagram number 458
    FFV1_0( w_fp[188], w_fp[169], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 459 OF 15495 ***
    // Wavefunction(s) for diagram number 459
    // (none)
    // Amplitude(s) for diagram number 459
    VVV1_0( w_fp[85], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];

    // *** DIAGRAM 460 OF 15495 ***
    // Wavefunction(s) for diagram number 460
    // (none)
    // Amplitude(s) for diagram number 460
    FFV1_0( w_fp[3], w_fp[190], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 461 OF 15495 ***
    // Wavefunction(s) for diagram number 461
    // (none)
    // Amplitude(s) for diagram number 461
    VVV1_0( w_fp[8], w_fp[194], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 462 OF 15495 ***
    // Wavefunction(s) for diagram number 462
    // (none)
    // Amplitude(s) for diagram number 462
    FFV1_0( w_fp[188], w_fp[190], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

    // *** DIAGRAM 463 OF 15495 ***
    // Wavefunction(s) for diagram number 463
    // (none)
    // Amplitude(s) for diagram number 463
    FFV1_0( w_fp[3], w_fp[169], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[92], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 464 OF 15495 ***
    // Wavefunction(s) for diagram number 464
    // (none)
    // Amplitude(s) for diagram number 464
    FFV1_0( w_fp[3], w_fp[197], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 465 OF 15495 ***
    // Wavefunction(s) for diagram number 465
    // (none)
    // Amplitude(s) for diagram number 465
    FFV1_0( w_fp[157], w_fp[169], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 466 OF 15495 ***
    // Wavefunction(s) for diagram number 466
    // (none)
    // Amplitude(s) for diagram number 466
    FFV1_0( w_fp[3], w_fp[169], w_fp[148], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[149], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[150], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 467 OF 15495 ***
    // Wavefunction(s) for diagram number 467
    FFV1_1( w_fp[2], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[197] );
    FFV1_1( w_fp[197], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[211] );
    // Amplitude(s) for diagram number 467
    FFV1_0( w_fp[162], w_fp[211], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 468 OF 15495 ***
    // Wavefunction(s) for diagram number 468
    // (none)
    // Amplitude(s) for diagram number 468
    FFV1_0( w_fp[160], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 469 OF 15495 ***
    // Wavefunction(s) for diagram number 469
    FFV1_1( w_fp[197], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[212] );
    // Amplitude(s) for diagram number 469
    FFV1_0( w_fp[192], w_fp[212], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 470 OF 15495 ***
    // Wavefunction(s) for diagram number 470
    // (none)
    // Amplitude(s) for diagram number 470
    FFV1_0( w_fp[160], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 471 OF 15495 ***
    // Wavefunction(s) for diagram number 471
    FFV1_1( w_fp[197], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[213] );
    // Amplitude(s) for diagram number 471
    FFV1_0( w_fp[192], w_fp[213], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 472 OF 15495 ***
    // Wavefunction(s) for diagram number 472
    // (none)
    // Amplitude(s) for diagram number 472
    FFV1_0( w_fp[162], w_fp[213], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 473 OF 15495 ***
    // Wavefunction(s) for diagram number 473
    FFV1P0_3( w_fp[3], w_fp[197], COUPs[1], 1.0, depCoup, 0., 0., w_fp[214] );
    // Amplitude(s) for diagram number 473
    VVVV1_0( w_fp[10], w_fp[214], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV3_0( w_fp[10], w_fp[214], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVVV4_0( w_fp[10], w_fp[214], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 474 OF 15495 ***
    // Wavefunction(s) for diagram number 474
    // (none)
    // Amplitude(s) for diagram number 474
    VVV1_0( w_fp[214], w_fp[7], w_fp[15], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 475 OF 15495 ***
    // Wavefunction(s) for diagram number 475
    // (none)
    // Amplitude(s) for diagram number 475
    VVV1_0( w_fp[214], w_fp[5], w_fp[13], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 476 OF 15495 ***
    // Wavefunction(s) for diagram number 476
    // (none)
    // Amplitude(s) for diagram number 476
    FFV1_0( w_fp[195], w_fp[212], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 477 OF 15495 ***
    // Wavefunction(s) for diagram number 477
    // (none)
    // Amplitude(s) for diagram number 477
    FFV1_0( w_fp[3], w_fp[212], w_fp[13], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 478 OF 15495 ***
    // Wavefunction(s) for diagram number 478
    // (none)
    // Amplitude(s) for diagram number 478
    FFV1_0( w_fp[195], w_fp[213], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 479 OF 15495 ***
    // Wavefunction(s) for diagram number 479
    // (none)
    // Amplitude(s) for diagram number 479
    FFV1_0( w_fp[3], w_fp[213], w_fp[15], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 480 OF 15495 ***
    // Wavefunction(s) for diagram number 480
    // (none)
    // Amplitude(s) for diagram number 480
    VVVV1_0( w_fp[26], w_fp[214], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    VVVV3_0( w_fp[26], w_fp[214], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    VVVV4_0( w_fp[26], w_fp[214], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];

    // *** DIAGRAM 481 OF 15495 ***
    // Wavefunction(s) for diagram number 481
    // (none)
    // Amplitude(s) for diagram number 481
    VVV1_0( w_fp[214], w_fp[7], w_fp[30], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];

    // *** DIAGRAM 482 OF 15495 ***
    // Wavefunction(s) for diagram number 482
    // (none)
    // Amplitude(s) for diagram number 482
    VVV1_0( w_fp[214], w_fp[4], w_fp[29], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];

    // *** DIAGRAM 483 OF 15495 ***
    // Wavefunction(s) for diagram number 483
    // (none)
    // Amplitude(s) for diagram number 483
    FFV1_0( w_fp[165], w_fp[211], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 484 OF 15495 ***
    // Wavefunction(s) for diagram number 484
    // (none)
    // Amplitude(s) for diagram number 484
    FFV1_0( w_fp[3], w_fp[211], w_fp[29], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 485 OF 15495 ***
    // Wavefunction(s) for diagram number 485
    // (none)
    // Amplitude(s) for diagram number 485
    FFV1_0( w_fp[165], w_fp[213], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];

    // *** DIAGRAM 486 OF 15495 ***
    // Wavefunction(s) for diagram number 486
    // (none)
    // Amplitude(s) for diagram number 486
    FFV1_0( w_fp[3], w_fp[213], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 487 OF 15495 ***
    // Wavefunction(s) for diagram number 487
    // (none)
    // Amplitude(s) for diagram number 487
    VVVV1_0( w_fp[44], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    VVVV3_0( w_fp[44], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    VVVV4_0( w_fp[44], w_fp[214], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 488 OF 15495 ***
    // Wavefunction(s) for diagram number 488
    // (none)
    // Amplitude(s) for diagram number 488
    VVV1_0( w_fp[214], w_fp[5], w_fp[47], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];

    // *** DIAGRAM 489 OF 15495 ***
    // Wavefunction(s) for diagram number 489
    // (none)
    // Amplitude(s) for diagram number 489
    VVV1_0( w_fp[214], w_fp[4], w_fp[45], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 490 OF 15495 ***
    // Wavefunction(s) for diagram number 490
    // (none)
    // Amplitude(s) for diagram number 490
    FFV1_0( w_fp[167], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];

    // *** DIAGRAM 491 OF 15495 ***
    // Wavefunction(s) for diagram number 491
    // (none)
    // Amplitude(s) for diagram number 491
    FFV1_0( w_fp[3], w_fp[211], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 492 OF 15495 ***
    // Wavefunction(s) for diagram number 492
    // (none)
    // Amplitude(s) for diagram number 492
    FFV1_0( w_fp[167], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 493 OF 15495 ***
    // Wavefunction(s) for diagram number 493
    // (none)
    // Amplitude(s) for diagram number 493
    FFV1_0( w_fp[3], w_fp[212], w_fp[47], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 494 OF 15495 ***
    // Wavefunction(s) for diagram number 494
    // (none)
    // Amplitude(s) for diagram number 494
    VVV1_0( w_fp[48], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    VVV1_0( w_fp[49], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    VVV1_0( w_fp[50], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 495 OF 15495 ***
    // Wavefunction(s) for diagram number 495
    // (none)
    // Amplitude(s) for diagram number 495
    FFV1_0( w_fp[3], w_fp[213], w_fp[48], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[213], w_fp[49], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[213], w_fp[50], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 496 OF 15495 ***
    // Wavefunction(s) for diagram number 496
    // (none)
    // Amplitude(s) for diagram number 496
    VVV1_0( w_fp[54], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    VVV1_0( w_fp[55], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    VVV1_0( w_fp[56], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 497 OF 15495 ***
    // Wavefunction(s) for diagram number 497
    // (none)
    // Amplitude(s) for diagram number 497
    FFV1_0( w_fp[3], w_fp[212], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[212], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[212], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 498 OF 15495 ***
    // Wavefunction(s) for diagram number 498
    // (none)
    // Amplitude(s) for diagram number 498
    VVV1_0( w_fp[60], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    VVV1_0( w_fp[61], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    VVV1_0( w_fp[62], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];

    // *** DIAGRAM 499 OF 15495 ***
    // Wavefunction(s) for diagram number 499
    // (none)
    // Amplitude(s) for diagram number 499
    FFV1_0( w_fp[3], w_fp[211], w_fp[60], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[211], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[211], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 500 OF 15495 ***
    // Wavefunction(s) for diagram number 500
    FFV1_1( w_fp[197], w_fp[8], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[215] );
    FFV1_2( w_fp[196], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[216] );
    // Amplitude(s) for diagram number 500
    FFV1_0( w_fp[216], w_fp[215], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 501 OF 15495 ***
    // Wavefunction(s) for diagram number 501
    // (none)
    // Amplitude(s) for diagram number 501
    FFV1_0( w_fp[199], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 502 OF 15495 ***
    // Wavefunction(s) for diagram number 502
    // (none)
    // Amplitude(s) for diagram number 502
    FFV1_0( w_fp[200], w_fp[212], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 503 OF 15495 ***
    // Wavefunction(s) for diagram number 503
    // (none)
    // Amplitude(s) for diagram number 503
    FFV1_0( w_fp[200], w_fp[213], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 504 OF 15495 ***
    // Wavefunction(s) for diagram number 504
    FFV1P0_3( w_fp[196], w_fp[197], COUPs[1], 1.0, depCoup, 0., 0., w_fp[217] );
    // Amplitude(s) for diagram number 504
    VVV1_0( w_fp[26], w_fp[217], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 505 OF 15495 ***
    // Wavefunction(s) for diagram number 505
    // (none)
    // Amplitude(s) for diagram number 505
    FFV1_0( w_fp[196], w_fp[213], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];

    // *** DIAGRAM 506 OF 15495 ***
    // Wavefunction(s) for diagram number 506
    // (none)
    // Amplitude(s) for diagram number 506
    FFV1_0( w_fp[199], w_fp[197], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];

    // *** DIAGRAM 507 OF 15495 ***
    // Wavefunction(s) for diagram number 507
    // (none)
    // Amplitude(s) for diagram number 507
    VVV1_0( w_fp[44], w_fp[217], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 508 OF 15495 ***
    // Wavefunction(s) for diagram number 508
    // (none)
    // Amplitude(s) for diagram number 508
    FFV1_0( w_fp[196], w_fp[212], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];

    // *** DIAGRAM 509 OF 15495 ***
    // Wavefunction(s) for diagram number 509
    // (none)
    // Amplitude(s) for diagram number 509
    FFV1_0( w_fp[216], w_fp[197], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 510 OF 15495 ***
    // Wavefunction(s) for diagram number 510
    // (none)
    // Amplitude(s) for diagram number 510
    FFV1_0( w_fp[199], w_fp[212], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 511 OF 15495 ***
    // Wavefunction(s) for diagram number 511
    // (none)
    // Amplitude(s) for diagram number 511
    FFV1_0( w_fp[216], w_fp[213], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 512 OF 15495 ***
    // Wavefunction(s) for diagram number 512
    // (none)
    // Amplitude(s) for diagram number 512
    FFV1_0( w_fp[196], w_fp[197], w_fp[60], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[197], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[197], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 513 OF 15495 ***
    // Wavefunction(s) for diagram number 513
    // (none)
    // Amplitude(s) for diagram number 513
    FFV1_0( w_fp[196], w_fp[215], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 514 OF 15495 ***
    // Wavefunction(s) for diagram number 514
    // (none)
    // Amplitude(s) for diagram number 514
    FFV1_0( w_fp[200], w_fp[197], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];

    // *** DIAGRAM 515 OF 15495 ***
    // Wavefunction(s) for diagram number 515
    // (none)
    // Amplitude(s) for diagram number 515
    FFV1_0( w_fp[196], w_fp[197], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 516 OF 15495 ***
    // Wavefunction(s) for diagram number 516
    FFV1_2( w_fp[168], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[218] );
    // Amplitude(s) for diagram number 516
    FFV1_0( w_fp[218], w_fp[215], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 517 OF 15495 ***
    // Wavefunction(s) for diagram number 517
    // (none)
    // Amplitude(s) for diagram number 517
    FFV1_0( w_fp[171], w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 518 OF 15495 ***
    // Wavefunction(s) for diagram number 518
    // (none)
    // Amplitude(s) for diagram number 518
    FFV1_0( w_fp[172], w_fp[211], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 519 OF 15495 ***
    // Wavefunction(s) for diagram number 519
    // (none)
    // Amplitude(s) for diagram number 519
    FFV1_0( w_fp[172], w_fp[213], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 520 OF 15495 ***
    // Wavefunction(s) for diagram number 520
    FFV1P0_3( w_fp[168], w_fp[197], COUPs[1], 1.0, depCoup, 0., 0., w_fp[219] );
    // Amplitude(s) for diagram number 520
    VVV1_0( w_fp[10], w_fp[219], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 521 OF 15495 ***
    // Wavefunction(s) for diagram number 521
    // (none)
    // Amplitude(s) for diagram number 521
    FFV1_0( w_fp[168], w_fp[213], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 522 OF 15495 ***
    // Wavefunction(s) for diagram number 522
    // (none)
    // Amplitude(s) for diagram number 522
    FFV1_0( w_fp[171], w_fp[197], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 523 OF 15495 ***
    // Wavefunction(s) for diagram number 523
    // (none)
    // Amplitude(s) for diagram number 523
    VVV1_0( w_fp[44], w_fp[219], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 524 OF 15495 ***
    // Wavefunction(s) for diagram number 524
    // (none)
    // Amplitude(s) for diagram number 524
    FFV1_0( w_fp[168], w_fp[211], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];

    // *** DIAGRAM 525 OF 15495 ***
    // Wavefunction(s) for diagram number 525
    // (none)
    // Amplitude(s) for diagram number 525
    FFV1_0( w_fp[218], w_fp[197], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];

    // *** DIAGRAM 526 OF 15495 ***
    // Wavefunction(s) for diagram number 526
    // (none)
    // Amplitude(s) for diagram number 526
    FFV1_0( w_fp[171], w_fp[211], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 527 OF 15495 ***
    // Wavefunction(s) for diagram number 527
    // (none)
    // Amplitude(s) for diagram number 527
    FFV1_0( w_fp[218], w_fp[213], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 528 OF 15495 ***
    // Wavefunction(s) for diagram number 528
    // (none)
    // Amplitude(s) for diagram number 528
    FFV1_0( w_fp[168], w_fp[197], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[197], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[197], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 529 OF 15495 ***
    // Wavefunction(s) for diagram number 529
    // (none)
    // Amplitude(s) for diagram number 529
    FFV1_0( w_fp[168], w_fp[215], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];

    // *** DIAGRAM 530 OF 15495 ***
    // Wavefunction(s) for diagram number 530
    // (none)
    // Amplitude(s) for diagram number 530
    FFV1_0( w_fp[172], w_fp[197], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 531 OF 15495 ***
    // Wavefunction(s) for diagram number 531
    // (none)
    // Amplitude(s) for diagram number 531
    FFV1_0( w_fp[168], w_fp[197], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 532 OF 15495 ***
    // Wavefunction(s) for diagram number 532
    // (none)
    // Amplitude(s) for diagram number 532
    FFV1_0( w_fp[204], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 533 OF 15495 ***
    // Wavefunction(s) for diagram number 533
    // (none)
    // Amplitude(s) for diagram number 533
    FFV1_0( w_fp[180], w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 534 OF 15495 ***
    // Wavefunction(s) for diagram number 534
    // (none)
    // Amplitude(s) for diagram number 534
    FFV1_0( w_fp[182], w_fp[211], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 535 OF 15495 ***
    // Wavefunction(s) for diagram number 535
    // (none)
    // Amplitude(s) for diagram number 535
    FFV1_0( w_fp[182], w_fp[212], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 536 OF 15495 ***
    // Wavefunction(s) for diagram number 536
    FFV1P0_3( w_fp[179], w_fp[197], COUPs[1], 1.0, depCoup, 0., 0., w_fp[220] );
    // Amplitude(s) for diagram number 536
    VVV1_0( w_fp[10], w_fp[220], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 537 OF 15495 ***
    // Wavefunction(s) for diagram number 537
    // (none)
    // Amplitude(s) for diagram number 537
    FFV1_0( w_fp[179], w_fp[212], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 538 OF 15495 ***
    // Wavefunction(s) for diagram number 538
    // (none)
    // Amplitude(s) for diagram number 538
    FFV1_0( w_fp[180], w_fp[197], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];

    // *** DIAGRAM 539 OF 15495 ***
    // Wavefunction(s) for diagram number 539
    // (none)
    // Amplitude(s) for diagram number 539
    VVV1_0( w_fp[26], w_fp[220], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 540 OF 15495 ***
    // Wavefunction(s) for diagram number 540
    // (none)
    // Amplitude(s) for diagram number 540
    FFV1_0( w_fp[179], w_fp[211], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];

    // *** DIAGRAM 541 OF 15495 ***
    // Wavefunction(s) for diagram number 541
    // (none)
    // Amplitude(s) for diagram number 541
    FFV1_0( w_fp[204], w_fp[197], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];

    // *** DIAGRAM 542 OF 15495 ***
    // Wavefunction(s) for diagram number 542
    // (none)
    // Amplitude(s) for diagram number 542
    FFV1_0( w_fp[180], w_fp[211], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 543 OF 15495 ***
    // Wavefunction(s) for diagram number 543
    // (none)
    // Amplitude(s) for diagram number 543
    FFV1_0( w_fp[204], w_fp[212], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 544 OF 15495 ***
    // Wavefunction(s) for diagram number 544
    // (none)
    // Amplitude(s) for diagram number 544
    FFV1_0( w_fp[179], w_fp[197], w_fp[48], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[197], w_fp[49], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[197], w_fp[50], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 545 OF 15495 ***
    // Wavefunction(s) for diagram number 545
    // (none)
    // Amplitude(s) for diagram number 545
    FFV1_0( w_fp[179], w_fp[215], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 546 OF 15495 ***
    // Wavefunction(s) for diagram number 546
    // (none)
    // Amplitude(s) for diagram number 546
    FFV1_0( w_fp[182], w_fp[197], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 547 OF 15495 ***
    // Wavefunction(s) for diagram number 547
    // (none)
    // Amplitude(s) for diagram number 547
    FFV1_0( w_fp[179], w_fp[197], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 548 OF 15495 ***
    // Wavefunction(s) for diagram number 548
    FFV1_2( w_fp[3], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[221] );
    // Amplitude(s) for diagram number 548
    FFV1_0( w_fp[221], w_fp[215], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 549 OF 15495 ***
    // Wavefunction(s) for diagram number 549
    // (none)
    // Amplitude(s) for diagram number 549
    FFV1_0( w_fp[3], w_fp[215], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 550 OF 15495 ***
    // Wavefunction(s) for diagram number 550
    FFV1_1( w_fp[197], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[222] );
    // Amplitude(s) for diagram number 550
    FFV1_0( w_fp[157], w_fp[222], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 551 OF 15495 ***
    // Wavefunction(s) for diagram number 551
    // (none)
    // Amplitude(s) for diagram number 551
    FFV1_0( w_fp[157], w_fp[213], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 552 OF 15495 ***
    // Wavefunction(s) for diagram number 552
    // (none)
    // Amplitude(s) for diagram number 552
    FFV1_0( w_fp[157], w_fp[197], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 553 OF 15495 ***
    // Wavefunction(s) for diagram number 553
    // (none)
    // Amplitude(s) for diagram number 553
    VVV1_0( w_fp[70], w_fp[214], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 554 OF 15495 ***
    // Wavefunction(s) for diagram number 554
    // (none)
    // Amplitude(s) for diagram number 554
    FFV1_0( w_fp[3], w_fp[213], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 555 OF 15495 ***
    // Wavefunction(s) for diagram number 555
    // (none)
    // Amplitude(s) for diagram number 555
    VVV1_0( w_fp[44], w_fp[214], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];

    // *** DIAGRAM 556 OF 15495 ***
    // Wavefunction(s) for diagram number 556
    // (none)
    // Amplitude(s) for diagram number 556
    FFV1_0( w_fp[3], w_fp[222], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 557 OF 15495 ***
    // Wavefunction(s) for diagram number 557
    // (none)
    // Amplitude(s) for diagram number 557
    FFV1_0( w_fp[221], w_fp[197], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 558 OF 15495 ***
    // Wavefunction(s) for diagram number 558
    // (none)
    // Amplitude(s) for diagram number 558
    VVV1_0( w_fp[8], w_fp[214], w_fp[69], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 559 OF 15495 ***
    // Wavefunction(s) for diagram number 559
    // (none)
    // Amplitude(s) for diagram number 559
    FFV1_0( w_fp[221], w_fp[213], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];

    // *** DIAGRAM 560 OF 15495 ***
    // Wavefunction(s) for diagram number 560
    // (none)
    // Amplitude(s) for diagram number 560
    FFV1_0( w_fp[3], w_fp[197], w_fp[81], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[82], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[83], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 561 OF 15495 ***
    // Wavefunction(s) for diagram number 561
    // (none)
    // Amplitude(s) for diagram number 561
    FFV1_0( w_fp[208], w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];

    // *** DIAGRAM 562 OF 15495 ***
    // Wavefunction(s) for diagram number 562
    // (none)
    // Amplitude(s) for diagram number 562
    FFV1_0( w_fp[3], w_fp[215], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 563 OF 15495 ***
    // Wavefunction(s) for diagram number 563
    FFV1_1( w_fp[197], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[223] );
    // Amplitude(s) for diagram number 563
    FFV1_0( w_fp[157], w_fp[223], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];

    // *** DIAGRAM 564 OF 15495 ***
    // Wavefunction(s) for diagram number 564
    // (none)
    // Amplitude(s) for diagram number 564
    FFV1_0( w_fp[157], w_fp[212], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 565 OF 15495 ***
    // Wavefunction(s) for diagram number 565
    // (none)
    // Amplitude(s) for diagram number 565
    FFV1_0( w_fp[157], w_fp[197], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 566 OF 15495 ***
    // Wavefunction(s) for diagram number 566
    // (none)
    // Amplitude(s) for diagram number 566
    VVV1_0( w_fp[105], w_fp[214], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 567 OF 15495 ***
    // Wavefunction(s) for diagram number 567
    // (none)
    // Amplitude(s) for diagram number 567
    FFV1_0( w_fp[3], w_fp[212], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 568 OF 15495 ***
    // Wavefunction(s) for diagram number 568
    // (none)
    // Amplitude(s) for diagram number 568
    VVV1_0( w_fp[26], w_fp[214], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];

    // *** DIAGRAM 569 OF 15495 ***
    // Wavefunction(s) for diagram number 569
    // (none)
    // Amplitude(s) for diagram number 569
    FFV1_0( w_fp[3], w_fp[223], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 570 OF 15495 ***
    // Wavefunction(s) for diagram number 570
    // (none)
    // Amplitude(s) for diagram number 570
    FFV1_0( w_fp[208], w_fp[197], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 571 OF 15495 ***
    // Wavefunction(s) for diagram number 571
    // (none)
    // Amplitude(s) for diagram number 571
    VVV1_0( w_fp[8], w_fp[214], w_fp[103], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];

    // *** DIAGRAM 572 OF 15495 ***
    // Wavefunction(s) for diagram number 572
    // (none)
    // Amplitude(s) for diagram number 572
    FFV1_0( w_fp[208], w_fp[212], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];

    // *** DIAGRAM 573 OF 15495 ***
    // Wavefunction(s) for diagram number 573
    // (none)
    // Amplitude(s) for diagram number 573
    FFV1_0( w_fp[3], w_fp[197], w_fp[107], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[108], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[109], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 574 OF 15495 ***
    // Wavefunction(s) for diagram number 574
    // (none)
    // Amplitude(s) for diagram number 574
    FFV1_0( w_fp[186], w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];

    // *** DIAGRAM 575 OF 15495 ***
    // Wavefunction(s) for diagram number 575
    // (none)
    // Amplitude(s) for diagram number 575
    FFV1_0( w_fp[3], w_fp[215], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 576 OF 15495 ***
    // Wavefunction(s) for diagram number 576
    // (none)
    // Amplitude(s) for diagram number 576
    FFV1_0( w_fp[157], w_fp[211], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 577 OF 15495 ***
    // Wavefunction(s) for diagram number 577
    FFV1_1( w_fp[197], w_fp[100], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[224] );
    // Amplitude(s) for diagram number 577
    FFV1_0( w_fp[157], w_fp[224], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 578 OF 15495 ***
    // Wavefunction(s) for diagram number 578
    // (none)
    // Amplitude(s) for diagram number 578
    FFV1_0( w_fp[157], w_fp[197], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 579 OF 15495 ***
    // Wavefunction(s) for diagram number 579
    // (none)
    // Amplitude(s) for diagram number 579
    VVV1_0( w_fp[10], w_fp[214], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 580 OF 15495 ***
    // Wavefunction(s) for diagram number 580
    // (none)
    // Amplitude(s) for diagram number 580
    FFV1_0( w_fp[3], w_fp[224], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 581 OF 15495 ***
    // Wavefunction(s) for diagram number 581
    // (none)
    // Amplitude(s) for diagram number 581
    FFV1_0( w_fp[186], w_fp[197], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 582 OF 15495 ***
    // Wavefunction(s) for diagram number 582
    // (none)
    // Amplitude(s) for diagram number 582
    VVV1_0( w_fp[101], w_fp[214], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];

    // *** DIAGRAM 583 OF 15495 ***
    // Wavefunction(s) for diagram number 583
    // (none)
    // Amplitude(s) for diagram number 583
    FFV1_0( w_fp[3], w_fp[211], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 584 OF 15495 ***
    // Wavefunction(s) for diagram number 584
    // (none)
    // Amplitude(s) for diagram number 584
    VVV1_0( w_fp[8], w_fp[214], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 585 OF 15495 ***
    // Wavefunction(s) for diagram number 585
    // (none)
    // Amplitude(s) for diagram number 585
    FFV1_0( w_fp[186], w_fp[211], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 586 OF 15495 ***
    // Wavefunction(s) for diagram number 586
    // (none)
    // Amplitude(s) for diagram number 586
    FFV1_0( w_fp[3], w_fp[197], w_fp[74], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[73], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[72], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 587 OF 15495 ***
    // Wavefunction(s) for diagram number 587
    // (none)
    // Amplitude(s) for diagram number 587
    FFV1_0( w_fp[3], w_fp[215], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 588 OF 15495 ***
    // Wavefunction(s) for diagram number 588
    // (none)
    // Amplitude(s) for diagram number 588
    FFV1_0( w_fp[157], w_fp[197], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[197], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 589 OF 15495 ***
    // Wavefunction(s) for diagram number 589
    // (none)
    // Amplitude(s) for diagram number 589
    FFV1_0( w_fp[3], w_fp[197], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[143], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[144], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 590 OF 15495 ***
    // Wavefunction(s) for diagram number 590
    FFV1_1( w_fp[2], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[215] );
    FFV1_1( w_fp[215], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[225] );
    // Amplitude(s) for diagram number 590
    FFV1_0( w_fp[162], w_fp[225], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 591 OF 15495 ***
    // Wavefunction(s) for diagram number 591
    // (none)
    // Amplitude(s) for diagram number 591
    FFV1_0( w_fp[159], w_fp[225], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 592 OF 15495 ***
    // Wavefunction(s) for diagram number 592
    FFV1_1( w_fp[215], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[226] );
    // Amplitude(s) for diagram number 592
    FFV1_0( w_fp[192], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 593 OF 15495 ***
    // Wavefunction(s) for diagram number 593
    // (none)
    // Amplitude(s) for diagram number 593
    FFV1_0( w_fp[159], w_fp[226], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 594 OF 15495 ***
    // Wavefunction(s) for diagram number 594
    FFV1_1( w_fp[215], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[227] );
    // Amplitude(s) for diagram number 594
    FFV1_0( w_fp[192], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 595 OF 15495 ***
    // Wavefunction(s) for diagram number 595
    // (none)
    // Amplitude(s) for diagram number 595
    FFV1_0( w_fp[162], w_fp[227], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 596 OF 15495 ***
    // Wavefunction(s) for diagram number 596
    FFV1P0_3( w_fp[3], w_fp[215], COUPs[1], 1.0, depCoup, 0., 0., w_fp[228] );
    // Amplitude(s) for diagram number 596
    VVVV1_0( w_fp[10], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[10], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV4_0( w_fp[10], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 597 OF 15495 ***
    // Wavefunction(s) for diagram number 597
    // (none)
    // Amplitude(s) for diagram number 597
    VVV1_0( w_fp[228], w_fp[6], w_fp[15], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 598 OF 15495 ***
    // Wavefunction(s) for diagram number 598
    // (none)
    // Amplitude(s) for diagram number 598
    VVV1_0( w_fp[228], w_fp[5], w_fp[12], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 599 OF 15495 ***
    // Wavefunction(s) for diagram number 599
    // (none)
    // Amplitude(s) for diagram number 599
    FFV1_0( w_fp[195], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 600 OF 15495 ***
    // Wavefunction(s) for diagram number 600
    // (none)
    // Amplitude(s) for diagram number 600
    FFV1_0( w_fp[3], w_fp[226], w_fp[12], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 197 );
    storeWf( wfs, w_cx, nevt, 204 );
    storeWf( wfs, w_cx, nevt, 205 );
    storeWf( wfs, w_cx, nevt, 206 );
    storeWf( wfs, w_cx, nevt, 207 );
    storeWf( wfs, w_cx, nevt, 208 );
    storeWf( wfs, w_cx, nevt, 209 );
    storeWf( wfs, w_cx, nevt, 210 );
    storeWf( wfs, w_cx, nevt, 211 );
    storeWf( wfs, w_cx, nevt, 212 );
    storeWf( wfs, w_cx, nevt, 213 );
    storeWf( wfs, w_cx, nevt, 214 );
    storeWf( wfs, w_cx, nevt, 215 );
    storeWf( wfs, w_cx, nevt, 216 );
    storeWf( wfs, w_cx, nevt, 217 );
    storeWf( wfs, w_cx, nevt, 218 );
    storeWf( wfs, w_cx, nevt, 219 );
    storeWf( wfs, w_cx, nevt, 220 );
    storeWf( wfs, w_cx, nevt, 221 );
    storeWf( wfs, w_cx, nevt, 222 );
    storeWf( wfs, w_cx, nevt, 223 );
    storeWf( wfs, w_cx, nevt, 224 );
    storeWf( wfs, w_cx, nevt, 225 );
    storeWf( wfs, w_cx, nevt, 226 );
    storeWf( wfs, w_cx, nevt, 227 );
    storeWf( wfs, w_cx, nevt, 228 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
