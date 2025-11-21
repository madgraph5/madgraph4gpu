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
  diagramgroup8( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 93 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 119 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 168 );
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
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 257 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 261 );
    retrieveWf( wfs, w_cx, nevt, 262 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 277 );
    retrieveWf( wfs, w_cx, nevt, 278 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 280 );
#endif
#endif

    // *** DIAGRAM 1401 OF 15495 ***
    // Wavefunction(s) for diagram number 1401
    // (none)
    // Amplitude(s) for diagram number 1401
    VVVV1_0( w_fp[258], w_fp[6], w_fp[7], w_fp[280], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[258], w_fp[6], w_fp[7], w_fp[280], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[258], w_fp[6], w_fp[7], w_fp[280], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1402 OF 15495 ***
    // Wavefunction(s) for diagram number 1402
    FFV1_2( w_fp[168], w_fp[258], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[281] );
    // Amplitude(s) for diagram number 1402
    FFV1_0( w_fp[281], w_fp[177], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1403 OF 15495 ***
    // Wavefunction(s) for diagram number 1403
    // (none)
    // Amplitude(s) for diagram number 1403
    FFV1_0( w_fp[168], w_fp[177], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];

    // *** DIAGRAM 1404 OF 15495 ***
    // Wavefunction(s) for diagram number 1404
    // (none)
    // Amplitude(s) for diagram number 1404
    FFV1_0( w_fp[171], w_fp[177], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1405 OF 15495 ***
    // Wavefunction(s) for diagram number 1405
    // (none)
    // Amplitude(s) for diagram number 1405
    FFV1_0( w_fp[281], w_fp[108], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1406 OF 15495 ***
    // Wavefunction(s) for diagram number 1406
    // (none)
    // Amplitude(s) for diagram number 1406
    FFV1_0( w_fp[168], w_fp[108], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1407 OF 15495 ***
    // Wavefunction(s) for diagram number 1407
    // (none)
    // Amplitude(s) for diagram number 1407
    FFV1_0( w_fp[170], w_fp[108], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1408 OF 15495 ***
    // Wavefunction(s) for diagram number 1408
    // (none)
    // Amplitude(s) for diagram number 1408
    FFV1_0( w_fp[171], w_fp[10], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];

    // *** DIAGRAM 1409 OF 15495 ***
    // Wavefunction(s) for diagram number 1409
    // (none)
    // Amplitude(s) for diagram number 1409
    FFV1_0( w_fp[170], w_fp[10], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 1410 OF 15495 ***
    // Wavefunction(s) for diagram number 1410
    // (none)
    // Amplitude(s) for diagram number 1410
    FFV1_0( w_fp[168], w_fp[279], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];

    // *** DIAGRAM 1411 OF 15495 ***
    // Wavefunction(s) for diagram number 1411
    // (none)
    // Amplitude(s) for diagram number 1411
    VVV1_0( w_fp[258], w_fp[84], w_fp[280], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1412 OF 15495 ***
    // Wavefunction(s) for diagram number 1412
    // (none)
    // Amplitude(s) for diagram number 1412
    FFV1_0( w_fp[168], w_fp[93], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1413 OF 15495 ***
    // Wavefunction(s) for diagram number 1413
    // (none)
    // Amplitude(s) for diagram number 1413
    FFV1_0( w_fp[175], w_fp[279], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1414 OF 15495 ***
    // Wavefunction(s) for diagram number 1414
    // (none)
    // Amplitude(s) for diagram number 1414
    FFV1_0( w_fp[176], w_fp[279], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1415 OF 15495 ***
    // Wavefunction(s) for diagram number 1415
    FFV1P0_3( w_fp[174], w_fp[10], COUPs[1], 1.0, depCoup, 0., 0., w_fp[282] );
    // Amplitude(s) for diagram number 1415
    VVV1_0( w_fp[282], w_fp[260], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1416 OF 15495 ***
    // Wavefunction(s) for diagram number 1416
    // (none)
    // Amplitude(s) for diagram number 1416
    VVV1_0( w_fp[282], w_fp[265], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1417 OF 15495 ***
    // Wavefunction(s) for diagram number 1417
    // (none)
    // Amplitude(s) for diagram number 1417
    VVVV1_0( w_fp[258], w_fp[5], w_fp[7], w_fp[282], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[258], w_fp[5], w_fp[7], w_fp[282], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[258], w_fp[5], w_fp[7], w_fp[282], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1418 OF 15495 ***
    // Wavefunction(s) for diagram number 1418
    FFV1_2( w_fp[174], w_fp[258], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[283] );
    // Amplitude(s) for diagram number 1418
    FFV1_0( w_fp[283], w_fp[95], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1419 OF 15495 ***
    // Wavefunction(s) for diagram number 1419
    // (none)
    // Amplitude(s) for diagram number 1419
    FFV1_0( w_fp[174], w_fp[95], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];

    // *** DIAGRAM 1420 OF 15495 ***
    // Wavefunction(s) for diagram number 1420
    // (none)
    // Amplitude(s) for diagram number 1420
    FFV1_0( w_fp[176], w_fp[95], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1421 OF 15495 ***
    // Wavefunction(s) for diagram number 1421
    // (none)
    // Amplitude(s) for diagram number 1421
    FFV1_0( w_fp[283], w_fp[108], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1422 OF 15495 ***
    // Wavefunction(s) for diagram number 1422
    // (none)
    // Amplitude(s) for diagram number 1422
    FFV1_0( w_fp[174], w_fp[108], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 1423 OF 15495 ***
    // Wavefunction(s) for diagram number 1423
    // (none)
    // Amplitude(s) for diagram number 1423
    FFV1_0( w_fp[175], w_fp[108], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1424 OF 15495 ***
    // Wavefunction(s) for diagram number 1424
    // (none)
    // Amplitude(s) for diagram number 1424
    FFV1_0( w_fp[176], w_fp[10], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];

    // *** DIAGRAM 1425 OF 15495 ***
    // Wavefunction(s) for diagram number 1425
    // (none)
    // Amplitude(s) for diagram number 1425
    FFV1_0( w_fp[175], w_fp[10], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];

    // *** DIAGRAM 1426 OF 15495 ***
    // Wavefunction(s) for diagram number 1426
    // (none)
    // Amplitude(s) for diagram number 1426
    FFV1_0( w_fp[174], w_fp[279], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];

    // *** DIAGRAM 1427 OF 15495 ***
    // Wavefunction(s) for diagram number 1427
    // (none)
    // Amplitude(s) for diagram number 1427
    VVV1_0( w_fp[258], w_fp[100], w_fp[282], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1428 OF 15495 ***
    // Wavefunction(s) for diagram number 1428
    // (none)
    // Amplitude(s) for diagram number 1428
    FFV1_0( w_fp[174], w_fp[97], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 1429 OF 15495 ***
    // Wavefunction(s) for diagram number 1429
    // (none)
    // Amplitude(s) for diagram number 1429
    FFV1_0( w_fp[180], w_fp[279], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1430 OF 15495 ***
    // Wavefunction(s) for diagram number 1430
    // (none)
    // Amplitude(s) for diagram number 1430
    FFV1_0( w_fp[181], w_fp[279], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1431 OF 15495 ***
    // Wavefunction(s) for diagram number 1431
    FFV1P0_3( w_fp[179], w_fp[10], COUPs[1], 1.0, depCoup, 0., 0., w_fp[284] );
    // Amplitude(s) for diagram number 1431
    VVV1_0( w_fp[284], w_fp[260], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1432 OF 15495 ***
    // Wavefunction(s) for diagram number 1432
    // (none)
    // Amplitude(s) for diagram number 1432
    VVV1_0( w_fp[284], w_fp[263], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1433 OF 15495 ***
    // Wavefunction(s) for diagram number 1433
    // (none)
    // Amplitude(s) for diagram number 1433
    VVVV1_0( w_fp[258], w_fp[5], w_fp[6], w_fp[284], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[258], w_fp[5], w_fp[6], w_fp[284], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[258], w_fp[5], w_fp[6], w_fp[284], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1434 OF 15495 ***
    // Wavefunction(s) for diagram number 1434
    FFV1_2( w_fp[179], w_fp[258], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[285] );
    // Amplitude(s) for diagram number 1434
    FFV1_0( w_fp[285], w_fp[95], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1435 OF 15495 ***
    // Wavefunction(s) for diagram number 1435
    // (none)
    // Amplitude(s) for diagram number 1435
    FFV1_0( w_fp[179], w_fp[95], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];

    // *** DIAGRAM 1436 OF 15495 ***
    // Wavefunction(s) for diagram number 1436
    // (none)
    // Amplitude(s) for diagram number 1436
    FFV1_0( w_fp[181], w_fp[95], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1437 OF 15495 ***
    // Wavefunction(s) for diagram number 1437
    // (none)
    // Amplitude(s) for diagram number 1437
    FFV1_0( w_fp[285], w_fp[177], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1438 OF 15495 ***
    // Wavefunction(s) for diagram number 1438
    // (none)
    // Amplitude(s) for diagram number 1438
    FFV1_0( w_fp[179], w_fp[177], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 1439 OF 15495 ***
    // Wavefunction(s) for diagram number 1439
    // (none)
    // Amplitude(s) for diagram number 1439
    FFV1_0( w_fp[180], w_fp[177], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1440 OF 15495 ***
    // Wavefunction(s) for diagram number 1440
    // (none)
    // Amplitude(s) for diagram number 1440
    FFV1_0( w_fp[181], w_fp[10], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];

    // *** DIAGRAM 1441 OF 15495 ***
    // Wavefunction(s) for diagram number 1441
    // (none)
    // Amplitude(s) for diagram number 1441
    FFV1_0( w_fp[180], w_fp[10], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];

    // *** DIAGRAM 1442 OF 15495 ***
    // Wavefunction(s) for diagram number 1442
    // (none)
    // Amplitude(s) for diagram number 1442
    FFV1_0( w_fp[179], w_fp[279], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];

    // *** DIAGRAM 1443 OF 15495 ***
    // Wavefunction(s) for diagram number 1443
    // (none)
    // Amplitude(s) for diagram number 1443
    VVV1_0( w_fp[258], w_fp[113], w_fp[284], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1444 OF 15495 ***
    // Wavefunction(s) for diagram number 1444
    // (none)
    // Amplitude(s) for diagram number 1444
    FFV1_0( w_fp[179], w_fp[127], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];

    // *** DIAGRAM 1445 OF 15495 ***
    // Wavefunction(s) for diagram number 1445
    // (none)
    // Amplitude(s) for diagram number 1445
    FFV1_0( w_fp[184], w_fp[279], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];

    // *** DIAGRAM 1446 OF 15495 ***
    // Wavefunction(s) for diagram number 1446
    // (none)
    // Amplitude(s) for diagram number 1446
    FFV1_0( w_fp[3], w_fp[279], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1447 OF 15495 ***
    // Wavefunction(s) for diagram number 1447
    VVV1P0_1( w_fp[258], w_fp[113], COUPs[0], 1.0, depCoup, 0., 0., w_fp[286] );
    // Amplitude(s) for diagram number 1447
    VVV1_0( w_fp[259], w_fp[286], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1448 OF 15495 ***
    // Wavefunction(s) for diagram number 1448
    // (none)
    // Amplitude(s) for diagram number 1448
    VVV1_0( w_fp[259], w_fp[265], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 1449 OF 15495 ***
    // Wavefunction(s) for diagram number 1449
    // (none)
    // Amplitude(s) for diagram number 1449
    VVV1_0( w_fp[259], w_fp[258], w_fp[116], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1450 OF 15495 ***
    // Wavefunction(s) for diagram number 1450
    // (none)
    // Amplitude(s) for diagram number 1450
    VVVV1_0( w_fp[258], w_fp[113], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV3_0( w_fp[258], w_fp[113], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    VVVV4_0( w_fp[258], w_fp[113], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 1451 OF 15495 ***
    // Wavefunction(s) for diagram number 1451
    // (none)
    // Amplitude(s) for diagram number 1451
    FFV1_0( w_fp[275], w_fp[127], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];

    // *** DIAGRAM 1452 OF 15495 ***
    // Wavefunction(s) for diagram number 1452
    // (none)
    // Amplitude(s) for diagram number 1452
    FFV1_0( w_fp[3], w_fp[127], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1453 OF 15495 ***
    // Wavefunction(s) for diagram number 1453
    // (none)
    // Amplitude(s) for diagram number 1453
    FFV1_0( w_fp[275], w_fp[108], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 1454 OF 15495 ***
    // Wavefunction(s) for diagram number 1454
    // (none)
    // Amplitude(s) for diagram number 1454
    FFV1_0( w_fp[3], w_fp[108], w_fp[286], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1455 OF 15495 ***
    // Wavefunction(s) for diagram number 1455
    // (none)
    // Amplitude(s) for diagram number 1455
    FFV1_0( w_fp[184], w_fp[108], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];

    // *** DIAGRAM 1456 OF 15495 ***
    // Wavefunction(s) for diagram number 1456
    // (none)
    // Amplitude(s) for diagram number 1456
    FFV1_0( w_fp[275], w_fp[10], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1457 OF 15495 ***
    // Wavefunction(s) for diagram number 1457
    // (none)
    // Amplitude(s) for diagram number 1457
    FFV1_0( w_fp[184], w_fp[10], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1458 OF 15495 ***
    // Wavefunction(s) for diagram number 1458
    // (none)
    // Amplitude(s) for diagram number 1458
    FFV1_0( w_fp[186], w_fp[279], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];

    // *** DIAGRAM 1459 OF 15495 ***
    // Wavefunction(s) for diagram number 1459
    // (none)
    // Amplitude(s) for diagram number 1459
    FFV1_0( w_fp[3], w_fp[279], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1460 OF 15495 ***
    // Wavefunction(s) for diagram number 1460
    VVV1P0_1( w_fp[258], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[287] );
    // Amplitude(s) for diagram number 1460
    VVV1_0( w_fp[259], w_fp[287], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];

    // *** DIAGRAM 1461 OF 15495 ***
    // Wavefunction(s) for diagram number 1461
    // (none)
    // Amplitude(s) for diagram number 1461
    VVV1_0( w_fp[259], w_fp[263], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 1462 OF 15495 ***
    // Wavefunction(s) for diagram number 1462
    // (none)
    // Amplitude(s) for diagram number 1462
    VVV1_0( w_fp[259], w_fp[258], w_fp[125], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 1463 OF 15495 ***
    // Wavefunction(s) for diagram number 1463
    // (none)
    // Amplitude(s) for diagram number 1463
    VVVV1_0( w_fp[258], w_fp[100], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    VVVV3_0( w_fp[258], w_fp[100], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    VVVV4_0( w_fp[258], w_fp[100], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];

    // *** DIAGRAM 1464 OF 15495 ***
    // Wavefunction(s) for diagram number 1464
    // (none)
    // Amplitude(s) for diagram number 1464
    FFV1_0( w_fp[275], w_fp[97], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];

    // *** DIAGRAM 1465 OF 15495 ***
    // Wavefunction(s) for diagram number 1465
    // (none)
    // Amplitude(s) for diagram number 1465
    FFV1_0( w_fp[3], w_fp[97], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1466 OF 15495 ***
    // Wavefunction(s) for diagram number 1466
    // (none)
    // Amplitude(s) for diagram number 1466
    FFV1_0( w_fp[275], w_fp[177], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];

    // *** DIAGRAM 1467 OF 15495 ***
    // Wavefunction(s) for diagram number 1467
    // (none)
    // Amplitude(s) for diagram number 1467
    FFV1_0( w_fp[3], w_fp[177], w_fp[287], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1468 OF 15495 ***
    // Wavefunction(s) for diagram number 1468
    // (none)
    // Amplitude(s) for diagram number 1468
    FFV1_0( w_fp[186], w_fp[177], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];

    // *** DIAGRAM 1469 OF 15495 ***
    // Wavefunction(s) for diagram number 1469
    // (none)
    // Amplitude(s) for diagram number 1469
    FFV1_0( w_fp[275], w_fp[10], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1470 OF 15495 ***
    // Wavefunction(s) for diagram number 1470
    // (none)
    // Amplitude(s) for diagram number 1470
    FFV1_0( w_fp[186], w_fp[10], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1471 OF 15495 ***
    // Wavefunction(s) for diagram number 1471
    // (none)
    // Amplitude(s) for diagram number 1471
    FFV1_0( w_fp[188], w_fp[279], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];

    // *** DIAGRAM 1472 OF 15495 ***
    // Wavefunction(s) for diagram number 1472
    // (none)
    // Amplitude(s) for diagram number 1472
    FFV1_0( w_fp[3], w_fp[279], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1473 OF 15495 ***
    // Wavefunction(s) for diagram number 1473
    // (none)
    // Amplitude(s) for diagram number 1473
    VVV1_0( w_fp[259], w_fp[260], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1474 OF 15495 ***
    // Wavefunction(s) for diagram number 1474
    VVV1P0_1( w_fp[258], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[288] );
    // Amplitude(s) for diagram number 1474
    VVV1_0( w_fp[259], w_fp[288], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];

    // *** DIAGRAM 1475 OF 15495 ***
    // Wavefunction(s) for diagram number 1475
    // (none)
    // Amplitude(s) for diagram number 1475
    VVV1_0( w_fp[259], w_fp[258], w_fp[131], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 1476 OF 15495 ***
    // Wavefunction(s) for diagram number 1476
    // (none)
    // Amplitude(s) for diagram number 1476
    VVVV1_0( w_fp[258], w_fp[5], w_fp[84], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV3_0( w_fp[258], w_fp[5], w_fp[84], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVVV4_0( w_fp[258], w_fp[5], w_fp[84], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 1477 OF 15495 ***
    // Wavefunction(s) for diagram number 1477
    // (none)
    // Amplitude(s) for diagram number 1477
    FFV1_0( w_fp[275], w_fp[95], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];

    // *** DIAGRAM 1478 OF 15495 ***
    // Wavefunction(s) for diagram number 1478
    // (none)
    // Amplitude(s) for diagram number 1478
    FFV1_0( w_fp[3], w_fp[95], w_fp[288], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1479 OF 15495 ***
    // Wavefunction(s) for diagram number 1479
    // (none)
    // Amplitude(s) for diagram number 1479
    FFV1_0( w_fp[188], w_fp[95], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];

    // *** DIAGRAM 1480 OF 15495 ***
    // Wavefunction(s) for diagram number 1480
    // (none)
    // Amplitude(s) for diagram number 1480
    FFV1_0( w_fp[275], w_fp[93], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 1481 OF 15495 ***
    // Wavefunction(s) for diagram number 1481
    // (none)
    // Amplitude(s) for diagram number 1481
    FFV1_0( w_fp[3], w_fp[93], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1482 OF 15495 ***
    // Wavefunction(s) for diagram number 1482
    // (none)
    // Amplitude(s) for diagram number 1482
    FFV1_0( w_fp[275], w_fp[10], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1483 OF 15495 ***
    // Wavefunction(s) for diagram number 1483
    // (none)
    // Amplitude(s) for diagram number 1483
    FFV1_0( w_fp[188], w_fp[10], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1484 OF 15495 ***
    // Wavefunction(s) for diagram number 1484
    // (none)
    // Amplitude(s) for diagram number 1484
    FFV1_0( w_fp[3], w_fp[279], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[279], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[279], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1485 OF 15495 ***
    // Wavefunction(s) for diagram number 1485
    // (none)
    // Amplitude(s) for diagram number 1485
    VVV1_0( w_fp[258], w_fp[151], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVV1_0( w_fp[258], w_fp[152], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    VVV1_0( w_fp[258], w_fp[153], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 1486 OF 15495 ***
    // Wavefunction(s) for diagram number 1486
    // (none)
    // Amplitude(s) for diagram number 1486
    FFV1_0( w_fp[3], w_fp[172], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[256], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[257], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1487 OF 15495 ***
    // Wavefunction(s) for diagram number 1487
    VVV1P0_1( w_fp[1], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[279] );
    VVV1P0_1( w_fp[279], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[289] );
    // Amplitude(s) for diagram number 1487
    VVVV1_0( w_fp[259], w_fp[289], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV3_0( w_fp[259], w_fp[289], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    VVVV4_0( w_fp[259], w_fp[289], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1488 OF 15495 ***
    // Wavefunction(s) for diagram number 1488
    // (none)
    // Amplitude(s) for diagram number 1488
    VVV1_0( w_fp[289], w_fp[7], w_fp[261], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];

    // *** DIAGRAM 1489 OF 15495 ***
    // Wavefunction(s) for diagram number 1489
    // (none)
    // Amplitude(s) for diagram number 1489
    VVV1_0( w_fp[289], w_fp[6], w_fp[262], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1490 OF 15495 ***
    // Wavefunction(s) for diagram number 1490
    VVV1P0_1( w_fp[279], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[290] );
    // Amplitude(s) for diagram number 1490
    VVVV1_0( w_fp[259], w_fp[290], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    VVVV3_0( w_fp[259], w_fp[290], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    VVVV4_0( w_fp[259], w_fp[290], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1491 OF 15495 ***
    // Wavefunction(s) for diagram number 1491
    VVV1P0_1( w_fp[259], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[291] );
    // Amplitude(s) for diagram number 1491
    VVV1_0( w_fp[290], w_fp[7], w_fp[291], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1492 OF 15495 ***
    // Wavefunction(s) for diagram number 1492
    // (none)
    // Amplitude(s) for diagram number 1492
    VVV1_0( w_fp[290], w_fp[4], w_fp[262], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1493 OF 15495 ***
    // Wavefunction(s) for diagram number 1493
    VVV1P0_1( w_fp[279], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[292] );
    // Amplitude(s) for diagram number 1493
    VVVV1_0( w_fp[259], w_fp[292], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    VVVV3_0( w_fp[259], w_fp[292], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    VVVV4_0( w_fp[259], w_fp[292], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];

    // *** DIAGRAM 1494 OF 15495 ***
    // Wavefunction(s) for diagram number 1494
    // (none)
    // Amplitude(s) for diagram number 1494
    VVV1_0( w_fp[292], w_fp[6], w_fp[291], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 1495 OF 15495 ***
    // Wavefunction(s) for diagram number 1495
    // (none)
    // Amplitude(s) for diagram number 1495
    VVV1_0( w_fp[292], w_fp[4], w_fp[261], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];

    // *** DIAGRAM 1496 OF 15495 ***
    // Wavefunction(s) for diagram number 1496
    VVVV1P0_1( w_fp[279], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[293] );
    VVVV3P0_1( w_fp[279], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[294] );
    VVVV4P0_1( w_fp[279], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[295] );
    // Amplitude(s) for diagram number 1496
    VVV1_0( w_fp[259], w_fp[293], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[294], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[295], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1497 OF 15495 ***
    // Wavefunction(s) for diagram number 1497
    VVVV1P0_1( w_fp[279], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[296] );
    VVVV3P0_1( w_fp[279], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[297] );
    VVVV4P0_1( w_fp[279], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[298] );
    // Amplitude(s) for diagram number 1497
    VVV1_0( w_fp[259], w_fp[296], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[297], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[298], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];

    // *** DIAGRAM 1498 OF 15495 ***
    // Wavefunction(s) for diagram number 1498
    VVVV1P0_1( w_fp[279], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[299] );
    VVVV3P0_1( w_fp[279], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[300] );
    VVVV4P0_1( w_fp[279], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[301] );
    // Amplitude(s) for diagram number 1498
    VVV1_0( w_fp[259], w_fp[299], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[300], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[301], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1499 OF 15495 ***
    // Wavefunction(s) for diagram number 1499
    FFV1_2( w_fp[3], w_fp[279], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[302] );
    // Amplitude(s) for diagram number 1499
    FFV1_0( w_fp[302], w_fp[63], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1500 OF 15495 ***
    // Wavefunction(s) for diagram number 1500
    // (none)
    // Amplitude(s) for diagram number 1500
    FFV1_0( w_fp[302], w_fp[46], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1501 OF 15495 ***
    // Wavefunction(s) for diagram number 1501
    FFV1P0_3( w_fp[3], w_fp[45], COUPs[1], 1.0, depCoup, 0., 0., w_fp[303] );
    // Amplitude(s) for diagram number 1501
    VVV1_0( w_fp[290], w_fp[7], w_fp[303], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1502 OF 15495 ***
    // Wavefunction(s) for diagram number 1502
    // (none)
    // Amplitude(s) for diagram number 1502
    FFV1_0( w_fp[3], w_fp[46], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];

    // *** DIAGRAM 1503 OF 15495 ***
    // Wavefunction(s) for diagram number 1503
    // (none)
    // Amplitude(s) for diagram number 1503
    VVV1_0( w_fp[292], w_fp[6], w_fp[303], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1504 OF 15495 ***
    // Wavefunction(s) for diagram number 1504
    // (none)
    // Amplitude(s) for diagram number 1504
    FFV1_0( w_fp[3], w_fp[63], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];

    // *** DIAGRAM 1505 OF 15495 ***
    // Wavefunction(s) for diagram number 1505
    // (none)
    // Amplitude(s) for diagram number 1505
    FFV1_0( w_fp[3], w_fp[45], w_fp[299], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[300], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[301], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1506 OF 15495 ***
    // Wavefunction(s) for diagram number 1506
    // (none)
    // Amplitude(s) for diagram number 1506
    FFV1_0( w_fp[302], w_fp[74], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1507 OF 15495 ***
    // Wavefunction(s) for diagram number 1507
    // (none)
    // Amplitude(s) for diagram number 1507
    FFV1_0( w_fp[302], w_fp[73], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1508 OF 15495 ***
    // Wavefunction(s) for diagram number 1508
    // (none)
    // Amplitude(s) for diagram number 1508
    VVV1_0( w_fp[289], w_fp[7], w_fp[277], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1509 OF 15495 ***
    // Wavefunction(s) for diagram number 1509
    // (none)
    // Amplitude(s) for diagram number 1509
    FFV1_0( w_fp[3], w_fp[73], w_fp[289], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];

    // *** DIAGRAM 1510 OF 15495 ***
    // Wavefunction(s) for diagram number 1510
    // (none)
    // Amplitude(s) for diagram number 1510
    VVV1_0( w_fp[292], w_fp[4], w_fp[277], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1511 OF 15495 ***
    // Wavefunction(s) for diagram number 1511
    // (none)
    // Amplitude(s) for diagram number 1511
    FFV1_0( w_fp[3], w_fp[74], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];

    // *** DIAGRAM 1512 OF 15495 ***
    // Wavefunction(s) for diagram number 1512
    // (none)
    // Amplitude(s) for diagram number 1512
    FFV1_0( w_fp[3], w_fp[177], w_fp[296], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[177], w_fp[297], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[177], w_fp[298], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1513 OF 15495 ***
    // Wavefunction(s) for diagram number 1513
    // (none)
    // Amplitude(s) for diagram number 1513
    FFV1_0( w_fp[302], w_fp[81], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1514 OF 15495 ***
    // Wavefunction(s) for diagram number 1514
    // (none)
    // Amplitude(s) for diagram number 1514
    FFV1_0( w_fp[302], w_fp[82], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1515 OF 15495 ***
    // Wavefunction(s) for diagram number 1515
    // (none)
    // Amplitude(s) for diagram number 1515
    VVV1_0( w_fp[289], w_fp[6], w_fp[278], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1516 OF 15495 ***
    // Wavefunction(s) for diagram number 1516
    // (none)
    // Amplitude(s) for diagram number 1516
    FFV1_0( w_fp[3], w_fp[82], w_fp[289], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 1517 OF 15495 ***
    // Wavefunction(s) for diagram number 1517
    // (none)
    // Amplitude(s) for diagram number 1517
    VVV1_0( w_fp[290], w_fp[4], w_fp[278], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1518 OF 15495 ***
    // Wavefunction(s) for diagram number 1518
    // (none)
    // Amplitude(s) for diagram number 1518
    FFV1_0( w_fp[3], w_fp[81], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 1519 OF 15495 ***
    // Wavefunction(s) for diagram number 1519
    // (none)
    // Amplitude(s) for diagram number 1519
    FFV1_0( w_fp[3], w_fp[108], w_fp[293], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[108], w_fp[294], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[108], w_fp[295], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1520 OF 15495 ***
    // Wavefunction(s) for diagram number 1520
    FFV1_1( w_fp[10], w_fp[279], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[304] );
    // Amplitude(s) for diagram number 1520
    FFV1_0( w_fp[198], w_fp[304], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1521 OF 15495 ***
    // Wavefunction(s) for diagram number 1521
    // (none)
    // Amplitude(s) for diagram number 1521
    FFV1_0( w_fp[199], w_fp[304], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1522 OF 15495 ***
    // Wavefunction(s) for diagram number 1522
    FFV1P0_3( w_fp[196], w_fp[10], COUPs[1], 1.0, depCoup, 0., 0., w_fp[305] );
    // Amplitude(s) for diagram number 1522
    VVV1_0( w_fp[305], w_fp[290], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1523 OF 15495 ***
    // Wavefunction(s) for diagram number 1523
    // (none)
    // Amplitude(s) for diagram number 1523
    VVV1_0( w_fp[305], w_fp[292], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1524 OF 15495 ***
    // Wavefunction(s) for diagram number 1524
    // (none)
    // Amplitude(s) for diagram number 1524
    VVVV1_0( w_fp[279], w_fp[6], w_fp[7], w_fp[305], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[279], w_fp[6], w_fp[7], w_fp[305], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[279], w_fp[6], w_fp[7], w_fp[305], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1525 OF 15495 ***
    // Wavefunction(s) for diagram number 1525
    FFV1_2( w_fp[196], w_fp[279], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[306] );
    // Amplitude(s) for diagram number 1525
    FFV1_0( w_fp[306], w_fp[177], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1526 OF 15495 ***
    // Wavefunction(s) for diagram number 1526
    // (none)
    // Amplitude(s) for diagram number 1526
    FFV1_0( w_fp[196], w_fp[177], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];

    // *** DIAGRAM 1527 OF 15495 ***
    // Wavefunction(s) for diagram number 1527
    // (none)
    // Amplitude(s) for diagram number 1527
    FFV1_0( w_fp[199], w_fp[177], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1528 OF 15495 ***
    // Wavefunction(s) for diagram number 1528
    // (none)
    // Amplitude(s) for diagram number 1528
    FFV1_0( w_fp[306], w_fp[108], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1529 OF 15495 ***
    // Wavefunction(s) for diagram number 1529
    // (none)
    // Amplitude(s) for diagram number 1529
    FFV1_0( w_fp[196], w_fp[108], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1530 OF 15495 ***
    // Wavefunction(s) for diagram number 1530
    // (none)
    // Amplitude(s) for diagram number 1530
    FFV1_0( w_fp[198], w_fp[108], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1531 OF 15495 ***
    // Wavefunction(s) for diagram number 1531
    // (none)
    // Amplitude(s) for diagram number 1531
    FFV1_0( w_fp[199], w_fp[10], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];

    // *** DIAGRAM 1532 OF 15495 ***
    // Wavefunction(s) for diagram number 1532
    // (none)
    // Amplitude(s) for diagram number 1532
    FFV1_0( w_fp[198], w_fp[10], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 1533 OF 15495 ***
    // Wavefunction(s) for diagram number 1533
    // (none)
    // Amplitude(s) for diagram number 1533
    FFV1_0( w_fp[196], w_fp[304], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];

    // *** DIAGRAM 1534 OF 15495 ***
    // Wavefunction(s) for diagram number 1534
    // (none)
    // Amplitude(s) for diagram number 1534
    VVV1_0( w_fp[279], w_fp[84], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1535 OF 15495 ***
    // Wavefunction(s) for diagram number 1535
    // (none)
    // Amplitude(s) for diagram number 1535
    FFV1_0( w_fp[196], w_fp[93], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1536 OF 15495 ***
    // Wavefunction(s) for diagram number 1536
    // (none)
    // Amplitude(s) for diagram number 1536
    FFV1_0( w_fp[202], w_fp[304], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1537 OF 15495 ***
    // Wavefunction(s) for diagram number 1537
    // (none)
    // Amplitude(s) for diagram number 1537
    FFV1_0( w_fp[176], w_fp[304], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1538 OF 15495 ***
    // Wavefunction(s) for diagram number 1538
    // (none)
    // Amplitude(s) for diagram number 1538
    VVV1_0( w_fp[282], w_fp[289], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1539 OF 15495 ***
    // Wavefunction(s) for diagram number 1539
    // (none)
    // Amplitude(s) for diagram number 1539
    VVV1_0( w_fp[282], w_fp[292], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1540 OF 15495 ***
    // Wavefunction(s) for diagram number 1540
    // (none)
    // Amplitude(s) for diagram number 1540
    VVVV1_0( w_fp[279], w_fp[4], w_fp[7], w_fp[282], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[279], w_fp[4], w_fp[7], w_fp[282], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[279], w_fp[4], w_fp[7], w_fp[282], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1541 OF 15495 ***
    // Wavefunction(s) for diagram number 1541
    FFV1_2( w_fp[174], w_fp[279], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[307] );
    // Amplitude(s) for diagram number 1541
    FFV1_0( w_fp[307], w_fp[45], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1542 OF 15495 ***
    // Wavefunction(s) for diagram number 1542
    // (none)
    // Amplitude(s) for diagram number 1542
    FFV1_0( w_fp[174], w_fp[45], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];

    // *** DIAGRAM 1543 OF 15495 ***
    // Wavefunction(s) for diagram number 1543
    // (none)
    // Amplitude(s) for diagram number 1543
    FFV1_0( w_fp[176], w_fp[45], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1544 OF 15495 ***
    // Wavefunction(s) for diagram number 1544
    // (none)
    // Amplitude(s) for diagram number 1544
    FFV1_0( w_fp[307], w_fp[108], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1545 OF 15495 ***
    // Wavefunction(s) for diagram number 1545
    // (none)
    // Amplitude(s) for diagram number 1545
    FFV1_0( w_fp[174], w_fp[108], w_fp[289], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];

    // *** DIAGRAM 1546 OF 15495 ***
    // Wavefunction(s) for diagram number 1546
    // (none)
    // Amplitude(s) for diagram number 1546
    FFV1_0( w_fp[202], w_fp[108], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1547 OF 15495 ***
    // Wavefunction(s) for diagram number 1547
    // (none)
    // Amplitude(s) for diagram number 1547
    FFV1_0( w_fp[176], w_fp[10], w_fp[289], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];

    // *** DIAGRAM 1548 OF 15495 ***
    // Wavefunction(s) for diagram number 1548
    // (none)
    // Amplitude(s) for diagram number 1548
    FFV1_0( w_fp[202], w_fp[10], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];

    // *** DIAGRAM 1549 OF 15495 ***
    // Wavefunction(s) for diagram number 1549
    // (none)
    // Amplitude(s) for diagram number 1549
    FFV1_0( w_fp[174], w_fp[304], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];

    // *** DIAGRAM 1550 OF 15495 ***
    // Wavefunction(s) for diagram number 1550
    // (none)
    // Amplitude(s) for diagram number 1550
    VVV1_0( w_fp[279], w_fp[102], w_fp[282], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1551 OF 15495 ***
    // Wavefunction(s) for diagram number 1551
    // (none)
    // Amplitude(s) for diagram number 1551
    FFV1_0( w_fp[174], w_fp[119], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];

    // *** DIAGRAM 1552 OF 15495 ***
    // Wavefunction(s) for diagram number 1552
    // (none)
    // Amplitude(s) for diagram number 1552
    FFV1_0( w_fp[204], w_fp[304], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1553 OF 15495 ***
    // Wavefunction(s) for diagram number 1553
    // (none)
    // Amplitude(s) for diagram number 1553
    FFV1_0( w_fp[181], w_fp[304], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1554 OF 15495 ***
    // Wavefunction(s) for diagram number 1554
    // (none)
    // Amplitude(s) for diagram number 1554
    VVV1_0( w_fp[284], w_fp[289], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1555 OF 15495 ***
    // Wavefunction(s) for diagram number 1555
    // (none)
    // Amplitude(s) for diagram number 1555
    VVV1_0( w_fp[284], w_fp[290], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1556 OF 15495 ***
    // Wavefunction(s) for diagram number 1556
    // (none)
    // Amplitude(s) for diagram number 1556
    VVVV1_0( w_fp[279], w_fp[4], w_fp[6], w_fp[284], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[279], w_fp[4], w_fp[6], w_fp[284], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[279], w_fp[4], w_fp[6], w_fp[284], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1557 OF 15495 ***
    // Wavefunction(s) for diagram number 1557
    FFV1_2( w_fp[179], w_fp[279], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[308] );
    // Amplitude(s) for diagram number 1557
    FFV1_0( w_fp[308], w_fp[45], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1558 OF 15495 ***
    // Wavefunction(s) for diagram number 1558
    // (none)
    // Amplitude(s) for diagram number 1558
    FFV1_0( w_fp[179], w_fp[45], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];

    // *** DIAGRAM 1559 OF 15495 ***
    // Wavefunction(s) for diagram number 1559
    // (none)
    // Amplitude(s) for diagram number 1559
    FFV1_0( w_fp[181], w_fp[45], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1560 OF 15495 ***
    // Wavefunction(s) for diagram number 1560
    // (none)
    // Amplitude(s) for diagram number 1560
    FFV1_0( w_fp[308], w_fp[177], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1561 OF 15495 ***
    // Wavefunction(s) for diagram number 1561
    // (none)
    // Amplitude(s) for diagram number 1561
    FFV1_0( w_fp[179], w_fp[177], w_fp[289], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];

    // *** DIAGRAM 1562 OF 15495 ***
    // Wavefunction(s) for diagram number 1562
    // (none)
    // Amplitude(s) for diagram number 1562
    FFV1_0( w_fp[204], w_fp[177], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1563 OF 15495 ***
    // Wavefunction(s) for diagram number 1563
    // (none)
    // Amplitude(s) for diagram number 1563
    FFV1_0( w_fp[181], w_fp[10], w_fp[289], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];

    // *** DIAGRAM 1564 OF 15495 ***
    // Wavefunction(s) for diagram number 1564
    // (none)
    // Amplitude(s) for diagram number 1564
    FFV1_0( w_fp[204], w_fp[10], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];

    // *** DIAGRAM 1565 OF 15495 ***
    // Wavefunction(s) for diagram number 1565
    // (none)
    // Amplitude(s) for diagram number 1565
    FFV1_0( w_fp[179], w_fp[304], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];

    // *** DIAGRAM 1566 OF 15495 ***
    // Wavefunction(s) for diagram number 1566
    // (none)
    // Amplitude(s) for diagram number 1566
    VVV1_0( w_fp[279], w_fp[86], w_fp[284], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1567 OF 15495 ***
    // Wavefunction(s) for diagram number 1567
    // (none)
    // Amplitude(s) for diagram number 1567
    FFV1_0( w_fp[179], w_fp[111], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];

    // *** DIAGRAM 1568 OF 15495 ***
    // Wavefunction(s) for diagram number 1568
    // (none)
    // Amplitude(s) for diagram number 1568
    FFV1_0( w_fp[206], w_fp[304], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];

    // *** DIAGRAM 1569 OF 15495 ***
    // Wavefunction(s) for diagram number 1569
    // (none)
    // Amplitude(s) for diagram number 1569
    FFV1_0( w_fp[3], w_fp[304], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1570 OF 15495 ***
    // Wavefunction(s) for diagram number 1570
    VVV1P0_1( w_fp[279], w_fp[86], COUPs[0], 1.0, depCoup, 0., 0., w_fp[309] );
    // Amplitude(s) for diagram number 1570
    VVV1_0( w_fp[259], w_fp[309], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 1571 OF 15495 ***
    // Wavefunction(s) for diagram number 1571
    // (none)
    // Amplitude(s) for diagram number 1571
    VVV1_0( w_fp[259], w_fp[292], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 1572 OF 15495 ***
    // Wavefunction(s) for diagram number 1572
    // (none)
    // Amplitude(s) for diagram number 1572
    VVV1_0( w_fp[259], w_fp[279], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 1573 OF 15495 ***
    // Wavefunction(s) for diagram number 1573
    // (none)
    // Amplitude(s) for diagram number 1573
    VVVV1_0( w_fp[279], w_fp[86], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVVV3_0( w_fp[279], w_fp[86], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    VVVV4_0( w_fp[279], w_fp[86], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];

    // *** DIAGRAM 1574 OF 15495 ***
    // Wavefunction(s) for diagram number 1574
    // (none)
    // Amplitude(s) for diagram number 1574
    FFV1_0( w_fp[302], w_fp[111], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];

    // *** DIAGRAM 1575 OF 15495 ***
    // Wavefunction(s) for diagram number 1575
    // (none)
    // Amplitude(s) for diagram number 1575
    FFV1_0( w_fp[3], w_fp[111], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1576 OF 15495 ***
    // Wavefunction(s) for diagram number 1576
    // (none)
    // Amplitude(s) for diagram number 1576
    FFV1_0( w_fp[302], w_fp[108], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];

    // *** DIAGRAM 1577 OF 15495 ***
    // Wavefunction(s) for diagram number 1577
    // (none)
    // Amplitude(s) for diagram number 1577
    FFV1_0( w_fp[3], w_fp[108], w_fp[309], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1578 OF 15495 ***
    // Wavefunction(s) for diagram number 1578
    // (none)
    // Amplitude(s) for diagram number 1578
    FFV1_0( w_fp[206], w_fp[108], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];

    // *** DIAGRAM 1579 OF 15495 ***
    // Wavefunction(s) for diagram number 1579
    // (none)
    // Amplitude(s) for diagram number 1579
    FFV1_0( w_fp[302], w_fp[10], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1580 OF 15495 ***
    // Wavefunction(s) for diagram number 1580
    // (none)
    // Amplitude(s) for diagram number 1580
    FFV1_0( w_fp[206], w_fp[10], w_fp[292], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1581 OF 15495 ***
    // Wavefunction(s) for diagram number 1581
    // (none)
    // Amplitude(s) for diagram number 1581
    FFV1_0( w_fp[208], w_fp[304], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];

    // *** DIAGRAM 1582 OF 15495 ***
    // Wavefunction(s) for diagram number 1582
    // (none)
    // Amplitude(s) for diagram number 1582
    FFV1_0( w_fp[3], w_fp[304], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1583 OF 15495 ***
    // Wavefunction(s) for diagram number 1583
    VVV1P0_1( w_fp[279], w_fp[102], COUPs[0], 1.0, depCoup, 0., 0., w_fp[310] );
    // Amplitude(s) for diagram number 1583
    VVV1_0( w_fp[259], w_fp[310], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];

    // *** DIAGRAM 1584 OF 15495 ***
    // Wavefunction(s) for diagram number 1584
    // (none)
    // Amplitude(s) for diagram number 1584
    VVV1_0( w_fp[259], w_fp[290], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 1585 OF 15495 ***
    // Wavefunction(s) for diagram number 1585
    // (none)
    // Amplitude(s) for diagram number 1585
    VVV1_0( w_fp[259], w_fp[279], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 1586 OF 15495 ***
    // Wavefunction(s) for diagram number 1586
    // (none)
    // Amplitude(s) for diagram number 1586
    VVVV1_0( w_fp[279], w_fp[102], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    VVVV3_0( w_fp[279], w_fp[102], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    VVVV4_0( w_fp[279], w_fp[102], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 1587 OF 15495 ***
    // Wavefunction(s) for diagram number 1587
    // (none)
    // Amplitude(s) for diagram number 1587
    FFV1_0( w_fp[302], w_fp[119], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 1588 OF 15495 ***
    // Wavefunction(s) for diagram number 1588
    // (none)
    // Amplitude(s) for diagram number 1588
    FFV1_0( w_fp[3], w_fp[119], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1589 OF 15495 ***
    // Wavefunction(s) for diagram number 1589
    // (none)
    // Amplitude(s) for diagram number 1589
    FFV1_0( w_fp[302], w_fp[177], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];

    // *** DIAGRAM 1590 OF 15495 ***
    // Wavefunction(s) for diagram number 1590
    // (none)
    // Amplitude(s) for diagram number 1590
    FFV1_0( w_fp[3], w_fp[177], w_fp[310], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1591 OF 15495 ***
    // Wavefunction(s) for diagram number 1591
    // (none)
    // Amplitude(s) for diagram number 1591
    FFV1_0( w_fp[208], w_fp[177], w_fp[279], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];

    // *** DIAGRAM 1592 OF 15495 ***
    // Wavefunction(s) for diagram number 1592
    // (none)
    // Amplitude(s) for diagram number 1592
    FFV1_0( w_fp[302], w_fp[10], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1593 OF 15495 ***
    // Wavefunction(s) for diagram number 1593
    // (none)
    // Amplitude(s) for diagram number 1593
    FFV1_0( w_fp[208], w_fp[10], w_fp[290], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1594 OF 15495 ***
    // Wavefunction(s) for diagram number 1594
    // (none)
    // Amplitude(s) for diagram number 1594
    FFV1_0( w_fp[188], w_fp[304], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];

    // *** DIAGRAM 1595 OF 15495 ***
    // Wavefunction(s) for diagram number 1595
    // (none)
    // Amplitude(s) for diagram number 1595
    FFV1_0( w_fp[3], w_fp[304], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1596 OF 15495 ***
    // Wavefunction(s) for diagram number 1596
    // (none)
    // Amplitude(s) for diagram number 1596
    VVV1_0( w_fp[259], w_fp[289], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1597 OF 15495 ***
    // Wavefunction(s) for diagram number 1597
    VVV1P0_1( w_fp[279], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[311] );
    // Amplitude(s) for diagram number 1597
    VVV1_0( w_fp[259], w_fp[311], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];

    // *** DIAGRAM 1598 OF 15495 ***
    // Wavefunction(s) for diagram number 1598
    // (none)
    // Amplitude(s) for diagram number 1598
    VVV1_0( w_fp[259], w_fp[279], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 1599 OF 15495 ***
    // Wavefunction(s) for diagram number 1599
    // (none)
    // Amplitude(s) for diagram number 1599
    VVVV1_0( w_fp[279], w_fp[4], w_fp[84], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV3_0( w_fp[279], w_fp[4], w_fp[84], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV4_0( w_fp[279], w_fp[4], w_fp[84], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];

    // *** DIAGRAM 1600 OF 15495 ***
    // Wavefunction(s) for diagram number 1600
    // (none)
    // Amplitude(s) for diagram number 1600
    FFV1_0( w_fp[302], w_fp[45], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 279 );
    storeWf( wfs, w_cx, nevt, 281 );
    storeWf( wfs, w_cx, nevt, 282 );
    storeWf( wfs, w_cx, nevt, 283 );
    storeWf( wfs, w_cx, nevt, 284 );
    storeWf( wfs, w_cx, nevt, 285 );
    storeWf( wfs, w_cx, nevt, 286 );
    storeWf( wfs, w_cx, nevt, 287 );
    storeWf( wfs, w_cx, nevt, 288 );
    storeWf( wfs, w_cx, nevt, 289 );
    storeWf( wfs, w_cx, nevt, 290 );
    storeWf( wfs, w_cx, nevt, 291 );
    storeWf( wfs, w_cx, nevt, 292 );
    storeWf( wfs, w_cx, nevt, 293 );
    storeWf( wfs, w_cx, nevt, 294 );
    storeWf( wfs, w_cx, nevt, 295 );
    storeWf( wfs, w_cx, nevt, 296 );
    storeWf( wfs, w_cx, nevt, 297 );
    storeWf( wfs, w_cx, nevt, 298 );
    storeWf( wfs, w_cx, nevt, 299 );
    storeWf( wfs, w_cx, nevt, 300 );
    storeWf( wfs, w_cx, nevt, 301 );
    storeWf( wfs, w_cx, nevt, 302 );
    storeWf( wfs, w_cx, nevt, 303 );
    storeWf( wfs, w_cx, nevt, 304 );
    storeWf( wfs, w_cx, nevt, 305 );
    storeWf( wfs, w_cx, nevt, 306 );
    storeWf( wfs, w_cx, nevt, 307 );
    storeWf( wfs, w_cx, nevt, 308 );
    storeWf( wfs, w_cx, nevt, 309 );
    storeWf( wfs, w_cx, nevt, 310 );
    storeWf( wfs, w_cx, nevt, 311 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
