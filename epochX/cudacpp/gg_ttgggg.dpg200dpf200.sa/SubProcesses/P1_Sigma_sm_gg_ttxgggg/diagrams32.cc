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
  diagramgroup32( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 156 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 158 );
    retrieveWf( wfs, w_cx, nevt, 163 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 457 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 466 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 517 );
    retrieveWf( wfs, w_cx, nevt, 534 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 549 );
#endif
#endif

    // *** DIAGRAM 6201 OF 15495 ***
    // Wavefunction(s) for diagram number 6201
    VVV1P0_1( w_fp[0], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[530] );
    FFV1_2( w_fp[3], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[532] );
    FFV1_2( w_fp[532], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[521] );
    // Amplitude(s) for diagram number 6201
    FFV1_0( w_fp[521], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6202 OF 15495 ***
    // Wavefunction(s) for diagram number 6202
    FFV1_2( w_fp[532], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[453] );
    // Amplitude(s) for diagram number 6202
    FFV1_0( w_fp[453], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6203 OF 15495 ***
    // Wavefunction(s) for diagram number 6203
    FFV1_2( w_fp[532], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[524] );
    // Amplitude(s) for diagram number 6203
    FFV1_0( w_fp[524], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6204 OF 15495 ***
    // Wavefunction(s) for diagram number 6204
    // (none)
    // Amplitude(s) for diagram number 6204
    FFV1_0( w_fp[453], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6205 OF 15495 ***
    // Wavefunction(s) for diagram number 6205
    // (none)
    // Amplitude(s) for diagram number 6205
    FFV1_0( w_fp[524], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6206 OF 15495 ***
    // Wavefunction(s) for diagram number 6206
    // (none)
    // Amplitude(s) for diagram number 6206
    FFV1_0( w_fp[521], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6207 OF 15495 ***
    // Wavefunction(s) for diagram number 6207
    VVV1P0_1( w_fp[530], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[523] );
    // Amplitude(s) for diagram number 6207
    VVVV1_0( w_fp[523], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[523], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    VVVV4_0( w_fp[523], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];

    // *** DIAGRAM 6208 OF 15495 ***
    // Wavefunction(s) for diagram number 6208
    VVV1P0_1( w_fp[523], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[580] );
    // Amplitude(s) for diagram number 6208
    VVV1_0( w_fp[456], w_fp[7], w_fp[580], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];

    // *** DIAGRAM 6209 OF 15495 ***
    // Wavefunction(s) for diagram number 6209
    VVV1P0_1( w_fp[523], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[515] );
    // Amplitude(s) for diagram number 6209
    VVV1_0( w_fp[456], w_fp[5], w_fp[515], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];

    // *** DIAGRAM 6210 OF 15495 ***
    // Wavefunction(s) for diagram number 6210
    FFV1_2( w_fp[3], w_fp[523], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[577] );
    // Amplitude(s) for diagram number 6210
    FFV1_0( w_fp[577], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];

    // *** DIAGRAM 6211 OF 15495 ***
    // Wavefunction(s) for diagram number 6211
    // (none)
    // Amplitude(s) for diagram number 6211
    FFV1_0( w_fp[3], w_fp[436], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6212 OF 15495 ***
    // Wavefunction(s) for diagram number 6212
    // (none)
    // Amplitude(s) for diagram number 6212
    FFV1_0( w_fp[577], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 6213 OF 15495 ***
    // Wavefunction(s) for diagram number 6213
    // (none)
    // Amplitude(s) for diagram number 6213
    FFV1_0( w_fp[3], w_fp[441], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6214 OF 15495 ***
    // Wavefunction(s) for diagram number 6214
    VVV1P0_1( w_fp[530], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[474] );
    // Amplitude(s) for diagram number 6214
    VVVV1_0( w_fp[474], w_fp[456], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[456], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[456], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6215 OF 15495 ***
    // Wavefunction(s) for diagram number 6215
    VVV1P0_1( w_fp[474], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[479] );
    // Amplitude(s) for diagram number 6215
    VVV1_0( w_fp[456], w_fp[7], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6216 OF 15495 ***
    // Wavefunction(s) for diagram number 6216
    VVV1P0_1( w_fp[474], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[536] );
    // Amplitude(s) for diagram number 6216
    VVV1_0( w_fp[456], w_fp[4], w_fp[536], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6217 OF 15495 ***
    // Wavefunction(s) for diagram number 6217
    FFV1_2( w_fp[3], w_fp[474], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[514] );
    // Amplitude(s) for diagram number 6217
    FFV1_0( w_fp[514], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];

    // *** DIAGRAM 6218 OF 15495 ***
    // Wavefunction(s) for diagram number 6218
    // (none)
    // Amplitude(s) for diagram number 6218
    FFV1_0( w_fp[3], w_fp[443], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6219 OF 15495 ***
    // Wavefunction(s) for diagram number 6219
    // (none)
    // Amplitude(s) for diagram number 6219
    FFV1_0( w_fp[514], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 6220 OF 15495 ***
    // Wavefunction(s) for diagram number 6220
    // (none)
    // Amplitude(s) for diagram number 6220
    FFV1_0( w_fp[3], w_fp[441], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6221 OF 15495 ***
    // Wavefunction(s) for diagram number 6221
    VVV1P0_1( w_fp[530], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[547] );
    // Amplitude(s) for diagram number 6221
    VVVV1_0( w_fp[547], w_fp[456], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    VVVV3_0( w_fp[547], w_fp[456], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    VVVV4_0( w_fp[547], w_fp[456], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6222 OF 15495 ***
    // Wavefunction(s) for diagram number 6222
    VVV1P0_1( w_fp[547], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[522] );
    // Amplitude(s) for diagram number 6222
    VVV1_0( w_fp[456], w_fp[5], w_fp[522], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];

    // *** DIAGRAM 6223 OF 15495 ***
    // Wavefunction(s) for diagram number 6223
    VVV1P0_1( w_fp[547], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[519] );
    // Amplitude(s) for diagram number 6223
    VVV1_0( w_fp[456], w_fp[4], w_fp[519], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6224 OF 15495 ***
    // Wavefunction(s) for diagram number 6224
    FFV1_2( w_fp[3], w_fp[547], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[488] );
    // Amplitude(s) for diagram number 6224
    FFV1_0( w_fp[488], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];

    // *** DIAGRAM 6225 OF 15495 ***
    // Wavefunction(s) for diagram number 6225
    // (none)
    // Amplitude(s) for diagram number 6225
    FFV1_0( w_fp[3], w_fp[443], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6226 OF 15495 ***
    // Wavefunction(s) for diagram number 6226
    // (none)
    // Amplitude(s) for diagram number 6226
    FFV1_0( w_fp[488], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];

    // *** DIAGRAM 6227 OF 15495 ***
    // Wavefunction(s) for diagram number 6227
    // (none)
    // Amplitude(s) for diagram number 6227
    FFV1_0( w_fp[3], w_fp[436], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6228 OF 15495 ***
    // Wavefunction(s) for diagram number 6228
    VVVV1P0_1( w_fp[530], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[588] );
    VVVV3P0_1( w_fp[530], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[586] );
    VVVV4P0_1( w_fp[530], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[451] );
    // Amplitude(s) for diagram number 6228
    VVV1_0( w_fp[588], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    VVV1_0( w_fp[586], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    VVV1_0( w_fp[451], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6229 OF 15495 ***
    // Wavefunction(s) for diagram number 6229
    // (none)
    // Amplitude(s) for diagram number 6229
    FFV1_0( w_fp[3], w_fp[441], w_fp[588], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[586], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6230 OF 15495 ***
    // Wavefunction(s) for diagram number 6230
    VVVV1P0_1( w_fp[530], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[520] );
    VVVV3P0_1( w_fp[530], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[486] );
    VVVV4P0_1( w_fp[530], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[576] );
    // Amplitude(s) for diagram number 6230
    VVV1_0( w_fp[520], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    VVV1_0( w_fp[486], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    VVV1_0( w_fp[576], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6231 OF 15495 ***
    // Wavefunction(s) for diagram number 6231
    // (none)
    // Amplitude(s) for diagram number 6231
    FFV1_0( w_fp[3], w_fp[436], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[436], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[436], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6232 OF 15495 ***
    // Wavefunction(s) for diagram number 6232
    VVVV1P0_1( w_fp[530], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[454] );
    VVVV3P0_1( w_fp[530], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[575] );
    VVVV4P0_1( w_fp[530], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[589] );
    // Amplitude(s) for diagram number 6232
    VVV1_0( w_fp[454], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    VVV1_0( w_fp[575], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    VVV1_0( w_fp[589], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 6233 OF 15495 ***
    // Wavefunction(s) for diagram number 6233
    // (none)
    // Amplitude(s) for diagram number 6233
    FFV1_0( w_fp[3], w_fp[443], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[443], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[443], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6234 OF 15495 ***
    // Wavefunction(s) for diagram number 6234
    FFV1_1( w_fp[259], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[442] );
    // Amplitude(s) for diagram number 6234
    FFV1_0( w_fp[216], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6235 OF 15495 ***
    // Wavefunction(s) for diagram number 6235
    // (none)
    // Amplitude(s) for diagram number 6235
    FFV1_0( w_fp[199], w_fp[442], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6236 OF 15495 ***
    // Wavefunction(s) for diagram number 6236
    FFV1_2( w_fp[196], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[529] );
    // Amplitude(s) for diagram number 6236
    FFV1_0( w_fp[529], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6237 OF 15495 ***
    // Wavefunction(s) for diagram number 6237
    // (none)
    // Amplitude(s) for diagram number 6237
    FFV1_0( w_fp[529], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6238 OF 15495 ***
    // Wavefunction(s) for diagram number 6238
    // (none)
    // Amplitude(s) for diagram number 6238
    VVV1_0( w_fp[474], w_fp[549], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6239 OF 15495 ***
    // Wavefunction(s) for diagram number 6239
    // (none)
    // Amplitude(s) for diagram number 6239
    FFV1_0( w_fp[196], w_fp[441], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 6240 OF 15495 ***
    // Wavefunction(s) for diagram number 6240
    // (none)
    // Amplitude(s) for diagram number 6240
    FFV1_0( w_fp[199], w_fp[259], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];

    // *** DIAGRAM 6241 OF 15495 ***
    // Wavefunction(s) for diagram number 6241
    // (none)
    // Amplitude(s) for diagram number 6241
    VVV1_0( w_fp[547], w_fp[549], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6242 OF 15495 ***
    // Wavefunction(s) for diagram number 6242
    // (none)
    // Amplitude(s) for diagram number 6242
    FFV1_0( w_fp[196], w_fp[436], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];

    // *** DIAGRAM 6243 OF 15495 ***
    // Wavefunction(s) for diagram number 6243
    // (none)
    // Amplitude(s) for diagram number 6243
    FFV1_0( w_fp[216], w_fp[259], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6244 OF 15495 ***
    // Wavefunction(s) for diagram number 6244
    // (none)
    // Amplitude(s) for diagram number 6244
    FFV1_0( w_fp[199], w_fp[436], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6245 OF 15495 ***
    // Wavefunction(s) for diagram number 6245
    // (none)
    // Amplitude(s) for diagram number 6245
    FFV1_0( w_fp[216], w_fp[441], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6246 OF 15495 ***
    // Wavefunction(s) for diagram number 6246
    // (none)
    // Amplitude(s) for diagram number 6246
    FFV1_0( w_fp[196], w_fp[259], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[259], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6247 OF 15495 ***
    // Wavefunction(s) for diagram number 6247
    // (none)
    // Amplitude(s) for diagram number 6247
    FFV1_0( w_fp[196], w_fp[442], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];

    // *** DIAGRAM 6248 OF 15495 ***
    // Wavefunction(s) for diagram number 6248
    // (none)
    // Amplitude(s) for diagram number 6248
    FFV1_0( w_fp[529], w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];

    // *** DIAGRAM 6249 OF 15495 ***
    // Wavefunction(s) for diagram number 6249
    VVV1P0_1( w_fp[530], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[510] );
    // Amplitude(s) for diagram number 6249
    FFV1_0( w_fp[196], w_fp[259], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6250 OF 15495 ***
    // Wavefunction(s) for diagram number 6250
    // (none)
    // Amplitude(s) for diagram number 6250
    FFV1_0( w_fp[218], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6251 OF 15495 ***
    // Wavefunction(s) for diagram number 6251
    // (none)
    // Amplitude(s) for diagram number 6251
    FFV1_0( w_fp[171], w_fp[442], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6252 OF 15495 ***
    // Wavefunction(s) for diagram number 6252
    FFV1_2( w_fp[168], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[555] );
    // Amplitude(s) for diagram number 6252
    FFV1_0( w_fp[555], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6253 OF 15495 ***
    // Wavefunction(s) for diagram number 6253
    // (none)
    // Amplitude(s) for diagram number 6253
    FFV1_0( w_fp[555], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6254 OF 15495 ***
    // Wavefunction(s) for diagram number 6254
    // (none)
    // Amplitude(s) for diagram number 6254
    VVV1_0( w_fp[523], w_fp[468], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6255 OF 15495 ***
    // Wavefunction(s) for diagram number 6255
    // (none)
    // Amplitude(s) for diagram number 6255
    FFV1_0( w_fp[168], w_fp[441], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];

    // *** DIAGRAM 6256 OF 15495 ***
    // Wavefunction(s) for diagram number 6256
    // (none)
    // Amplitude(s) for diagram number 6256
    FFV1_0( w_fp[171], w_fp[259], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];

    // *** DIAGRAM 6257 OF 15495 ***
    // Wavefunction(s) for diagram number 6257
    // (none)
    // Amplitude(s) for diagram number 6257
    VVV1_0( w_fp[547], w_fp[468], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6258 OF 15495 ***
    // Wavefunction(s) for diagram number 6258
    // (none)
    // Amplitude(s) for diagram number 6258
    FFV1_0( w_fp[168], w_fp[443], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];

    // *** DIAGRAM 6259 OF 15495 ***
    // Wavefunction(s) for diagram number 6259
    // (none)
    // Amplitude(s) for diagram number 6259
    FFV1_0( w_fp[218], w_fp[259], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];

    // *** DIAGRAM 6260 OF 15495 ***
    // Wavefunction(s) for diagram number 6260
    // (none)
    // Amplitude(s) for diagram number 6260
    FFV1_0( w_fp[171], w_fp[443], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6261 OF 15495 ***
    // Wavefunction(s) for diagram number 6261
    // (none)
    // Amplitude(s) for diagram number 6261
    FFV1_0( w_fp[218], w_fp[441], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6262 OF 15495 ***
    // Wavefunction(s) for diagram number 6262
    // (none)
    // Amplitude(s) for diagram number 6262
    FFV1_0( w_fp[168], w_fp[259], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6263 OF 15495 ***
    // Wavefunction(s) for diagram number 6263
    // (none)
    // Amplitude(s) for diagram number 6263
    FFV1_0( w_fp[168], w_fp[442], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];

    // *** DIAGRAM 6264 OF 15495 ***
    // Wavefunction(s) for diagram number 6264
    // (none)
    // Amplitude(s) for diagram number 6264
    FFV1_0( w_fp[555], w_fp[259], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];

    // *** DIAGRAM 6265 OF 15495 ***
    // Wavefunction(s) for diagram number 6265
    VVV1P0_1( w_fp[530], w_fp[102], COUPs[0], 1.0, depCoup, 0., 0., w_fp[556] );
    // Amplitude(s) for diagram number 6265
    FFV1_0( w_fp[168], w_fp[259], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6266 OF 15495 ***
    // Wavefunction(s) for diagram number 6266
    // (none)
    // Amplitude(s) for diagram number 6266
    FFV1_0( w_fp[204], w_fp[442], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6267 OF 15495 ***
    // Wavefunction(s) for diagram number 6267
    // (none)
    // Amplitude(s) for diagram number 6267
    FFV1_0( w_fp[180], w_fp[442], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6268 OF 15495 ***
    // Wavefunction(s) for diagram number 6268
    FFV1_2( w_fp[179], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[557] );
    // Amplitude(s) for diagram number 6268
    FFV1_0( w_fp[557], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6269 OF 15495 ***
    // Wavefunction(s) for diagram number 6269
    // (none)
    // Amplitude(s) for diagram number 6269
    FFV1_0( w_fp[557], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6270 OF 15495 ***
    // Wavefunction(s) for diagram number 6270
    // (none)
    // Amplitude(s) for diagram number 6270
    VVV1_0( w_fp[523], w_fp[517], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6271 OF 15495 ***
    // Wavefunction(s) for diagram number 6271
    // (none)
    // Amplitude(s) for diagram number 6271
    FFV1_0( w_fp[179], w_fp[436], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];

    // *** DIAGRAM 6272 OF 15495 ***
    // Wavefunction(s) for diagram number 6272
    // (none)
    // Amplitude(s) for diagram number 6272
    FFV1_0( w_fp[180], w_fp[259], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];

    // *** DIAGRAM 6273 OF 15495 ***
    // Wavefunction(s) for diagram number 6273
    // (none)
    // Amplitude(s) for diagram number 6273
    VVV1_0( w_fp[474], w_fp[517], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6274 OF 15495 ***
    // Wavefunction(s) for diagram number 6274
    // (none)
    // Amplitude(s) for diagram number 6274
    FFV1_0( w_fp[179], w_fp[443], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];

    // *** DIAGRAM 6275 OF 15495 ***
    // Wavefunction(s) for diagram number 6275
    // (none)
    // Amplitude(s) for diagram number 6275
    FFV1_0( w_fp[204], w_fp[259], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];

    // *** DIAGRAM 6276 OF 15495 ***
    // Wavefunction(s) for diagram number 6276
    // (none)
    // Amplitude(s) for diagram number 6276
    FFV1_0( w_fp[180], w_fp[443], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6277 OF 15495 ***
    // Wavefunction(s) for diagram number 6277
    // (none)
    // Amplitude(s) for diagram number 6277
    FFV1_0( w_fp[204], w_fp[436], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6278 OF 15495 ***
    // Wavefunction(s) for diagram number 6278
    // (none)
    // Amplitude(s) for diagram number 6278
    FFV1_0( w_fp[179], w_fp[259], w_fp[588], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[586], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6279 OF 15495 ***
    // Wavefunction(s) for diagram number 6279
    // (none)
    // Amplitude(s) for diagram number 6279
    FFV1_0( w_fp[179], w_fp[442], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];

    // *** DIAGRAM 6280 OF 15495 ***
    // Wavefunction(s) for diagram number 6280
    // (none)
    // Amplitude(s) for diagram number 6280
    FFV1_0( w_fp[557], w_fp[259], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];

    // *** DIAGRAM 6281 OF 15495 ***
    // Wavefunction(s) for diagram number 6281
    VVV1P0_1( w_fp[530], w_fp[66], COUPs[0], 1.0, depCoup, 0., 0., w_fp[518] );
    // Amplitude(s) for diagram number 6281
    FFV1_0( w_fp[179], w_fp[259], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6282 OF 15495 ***
    // Wavefunction(s) for diagram number 6282
    // (none)
    // Amplitude(s) for diagram number 6282
    FFV1_0( w_fp[221], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];

    // *** DIAGRAM 6283 OF 15495 ***
    // Wavefunction(s) for diagram number 6283
    // (none)
    // Amplitude(s) for diagram number 6283
    FFV1_0( w_fp[3], w_fp[442], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6284 OF 15495 ***
    // Wavefunction(s) for diagram number 6284
    // (none)
    // Amplitude(s) for diagram number 6284
    FFV1_0( w_fp[532], w_fp[457], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];

    // *** DIAGRAM 6285 OF 15495 ***
    // Wavefunction(s) for diagram number 6285
    // (none)
    // Amplitude(s) for diagram number 6285
    FFV1_0( w_fp[532], w_fp[441], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6286 OF 15495 ***
    // Wavefunction(s) for diagram number 6286
    // (none)
    // Amplitude(s) for diagram number 6286
    FFV1_0( w_fp[532], w_fp[259], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6287 OF 15495 ***
    // Wavefunction(s) for diagram number 6287
    // (none)
    // Amplitude(s) for diagram number 6287
    VVV1_0( w_fp[518], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 6288 OF 15495 ***
    // Wavefunction(s) for diagram number 6288
    // (none)
    // Amplitude(s) for diagram number 6288
    FFV1_0( w_fp[3], w_fp[441], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6289 OF 15495 ***
    // Wavefunction(s) for diagram number 6289
    // (none)
    // Amplitude(s) for diagram number 6289
    VVV1_0( w_fp[547], w_fp[456], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 6290 OF 15495 ***
    // Wavefunction(s) for diagram number 6290
    // (none)
    // Amplitude(s) for diagram number 6290
    FFV1_0( w_fp[3], w_fp[457], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6291 OF 15495 ***
    // Wavefunction(s) for diagram number 6291
    // (none)
    // Amplitude(s) for diagram number 6291
    FFV1_0( w_fp[221], w_fp[259], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6292 OF 15495 ***
    // Wavefunction(s) for diagram number 6292
    // (none)
    // Amplitude(s) for diagram number 6292
    VVV1_0( w_fp[530], w_fp[456], w_fp[69], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6293 OF 15495 ***
    // Wavefunction(s) for diagram number 6293
    // (none)
    // Amplitude(s) for diagram number 6293
    FFV1_0( w_fp[221], w_fp[441], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 6294 OF 15495 ***
    // Wavefunction(s) for diagram number 6294
    VVVV1P0_1( w_fp[530], w_fp[66], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[585] );
    VVVV3P0_1( w_fp[530], w_fp[66], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[541] );
    VVVV4P0_1( w_fp[530], w_fp[66], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[546] );
    // Amplitude(s) for diagram number 6294
    FFV1_0( w_fp[3], w_fp[259], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[541], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 6295 OF 15495 ***
    // Wavefunction(s) for diagram number 6295
    // (none)
    // Amplitude(s) for diagram number 6295
    FFV1_0( w_fp[208], w_fp[442], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];

    // *** DIAGRAM 6296 OF 15495 ***
    // Wavefunction(s) for diagram number 6296
    // (none)
    // Amplitude(s) for diagram number 6296
    FFV1_0( w_fp[3], w_fp[442], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6297 OF 15495 ***
    // Wavefunction(s) for diagram number 6297
    // (none)
    // Amplitude(s) for diagram number 6297
    FFV1_0( w_fp[532], w_fp[463], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];

    // *** DIAGRAM 6298 OF 15495 ***
    // Wavefunction(s) for diagram number 6298
    // (none)
    // Amplitude(s) for diagram number 6298
    FFV1_0( w_fp[532], w_fp[436], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];

    // *** DIAGRAM 6299 OF 15495 ***
    // Wavefunction(s) for diagram number 6299
    // (none)
    // Amplitude(s) for diagram number 6299
    FFV1_0( w_fp[532], w_fp[259], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6300 OF 15495 ***
    // Wavefunction(s) for diagram number 6300
    // (none)
    // Amplitude(s) for diagram number 6300
    VVV1_0( w_fp[556], w_fp[456], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 6301 OF 15495 ***
    // Wavefunction(s) for diagram number 6301
    // (none)
    // Amplitude(s) for diagram number 6301
    FFV1_0( w_fp[3], w_fp[436], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6302 OF 15495 ***
    // Wavefunction(s) for diagram number 6302
    // (none)
    // Amplitude(s) for diagram number 6302
    VVV1_0( w_fp[474], w_fp[456], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 6303 OF 15495 ***
    // Wavefunction(s) for diagram number 6303
    // (none)
    // Amplitude(s) for diagram number 6303
    FFV1_0( w_fp[3], w_fp[463], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6304 OF 15495 ***
    // Wavefunction(s) for diagram number 6304
    // (none)
    // Amplitude(s) for diagram number 6304
    FFV1_0( w_fp[208], w_fp[259], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6305 OF 15495 ***
    // Wavefunction(s) for diagram number 6305
    // (none)
    // Amplitude(s) for diagram number 6305
    VVV1_0( w_fp[530], w_fp[456], w_fp[103], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];

    // *** DIAGRAM 6306 OF 15495 ***
    // Wavefunction(s) for diagram number 6306
    // (none)
    // Amplitude(s) for diagram number 6306
    FFV1_0( w_fp[208], w_fp[436], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];

    // *** DIAGRAM 6307 OF 15495 ***
    // Wavefunction(s) for diagram number 6307
    VVVV1P0_1( w_fp[530], w_fp[102], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[550] );
    VVVV3P0_1( w_fp[530], w_fp[102], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[480] );
    VVVV4P0_1( w_fp[530], w_fp[102], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[558] );
    // Amplitude(s) for diagram number 6307
    FFV1_0( w_fp[3], w_fp[259], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[226] -= amp_sv[0];

    // *** DIAGRAM 6308 OF 15495 ***
    // Wavefunction(s) for diagram number 6308
    // (none)
    // Amplitude(s) for diagram number 6308
    FFV1_0( w_fp[186], w_fp[442], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];

    // *** DIAGRAM 6309 OF 15495 ***
    // Wavefunction(s) for diagram number 6309
    // (none)
    // Amplitude(s) for diagram number 6309
    FFV1_0( w_fp[3], w_fp[442], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6310 OF 15495 ***
    // Wavefunction(s) for diagram number 6310
    // (none)
    // Amplitude(s) for diagram number 6310
    FFV1_0( w_fp[532], w_fp[443], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];

    // *** DIAGRAM 6311 OF 15495 ***
    // Wavefunction(s) for diagram number 6311
    // (none)
    // Amplitude(s) for diagram number 6311
    FFV1_0( w_fp[532], w_fp[466], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6312 OF 15495 ***
    // Wavefunction(s) for diagram number 6312
    // (none)
    // Amplitude(s) for diagram number 6312
    FFV1_0( w_fp[532], w_fp[259], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6313 OF 15495 ***
    // Wavefunction(s) for diagram number 6313
    // (none)
    // Amplitude(s) for diagram number 6313
    VVV1_0( w_fp[523], w_fp[456], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6314 OF 15495 ***
    // Wavefunction(s) for diagram number 6314
    // (none)
    // Amplitude(s) for diagram number 6314
    FFV1_0( w_fp[3], w_fp[466], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6315 OF 15495 ***
    // Wavefunction(s) for diagram number 6315
    // (none)
    // Amplitude(s) for diagram number 6315
    FFV1_0( w_fp[186], w_fp[259], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6316 OF 15495 ***
    // Wavefunction(s) for diagram number 6316
    // (none)
    // Amplitude(s) for diagram number 6316
    VVV1_0( w_fp[510], w_fp[456], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];

    // *** DIAGRAM 6317 OF 15495 ***
    // Wavefunction(s) for diagram number 6317
    // (none)
    // Amplitude(s) for diagram number 6317
    FFV1_0( w_fp[3], w_fp[443], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6318 OF 15495 ***
    // Wavefunction(s) for diagram number 6318
    // (none)
    // Amplitude(s) for diagram number 6318
    VVV1_0( w_fp[530], w_fp[456], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6319 OF 15495 ***
    // Wavefunction(s) for diagram number 6319
    // (none)
    // Amplitude(s) for diagram number 6319
    FFV1_0( w_fp[186], w_fp[443], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];

    // *** DIAGRAM 6320 OF 15495 ***
    // Wavefunction(s) for diagram number 6320
    VVVV1P0_1( w_fp[530], w_fp[4], w_fp[100], COUPs[2], 1.0, depCoup, 0., 0., w_fp[559] );
    VVVV3P0_1( w_fp[530], w_fp[4], w_fp[100], COUPs[2], 1.0, depCoup, 0., 0., w_fp[560] );
    VVVV4P0_1( w_fp[530], w_fp[4], w_fp[100], COUPs[2], 1.0, depCoup, 0., 0., w_fp[438] );
    // Amplitude(s) for diagram number 6320
    FFV1_0( w_fp[3], w_fp[259], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6321 OF 15495 ***
    // Wavefunction(s) for diagram number 6321
    // (none)
    // Amplitude(s) for diagram number 6321
    FFV1_0( w_fp[3], w_fp[442], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[442], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[442], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6322 OF 15495 ***
    // Wavefunction(s) for diagram number 6322
    // (none)
    // Amplitude(s) for diagram number 6322
    FFV1_0( w_fp[532], w_fp[259], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[532], w_fp[259], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[532], w_fp[259], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6323 OF 15495 ***
    // Wavefunction(s) for diagram number 6323
    VVV1P0_1( w_fp[530], w_fp[139], COUPs[0], 1.0, depCoup, 0., 0., w_fp[442] );
    VVV1P0_1( w_fp[530], w_fp[140], COUPs[0], 1.0, depCoup, 0., 0., w_fp[562] );
    VVV1P0_1( w_fp[530], w_fp[141], COUPs[0], 1.0, depCoup, 0., 0., w_fp[437] );
    // Amplitude(s) for diagram number 6323
    FFV1_0( w_fp[3], w_fp[259], w_fp[442], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[133] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[165] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[437], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 6324 OF 15495 ***
    // Wavefunction(s) for diagram number 6324
    FFV1_1( w_fp[2], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[449] );
    FFV1_1( w_fp[449], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[590] );
    // Amplitude(s) for diagram number 6324
    FFV1_0( w_fp[94], w_fp[590], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6325 OF 15495 ***
    // Wavefunction(s) for diagram number 6325
    FFV1_1( w_fp[449], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[481] );
    // Amplitude(s) for diagram number 6325
    FFV1_0( w_fp[94], w_fp[481], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6326 OF 15495 ***
    // Wavefunction(s) for diagram number 6326
    FFV1_1( w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[596] );
    // Amplitude(s) for diagram number 6326
    FFV1_0( w_fp[29], w_fp[596], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6327 OF 15495 ***
    // Wavefunction(s) for diagram number 6327
    // (none)
    // Amplitude(s) for diagram number 6327
    FFV1_0( w_fp[29], w_fp[481], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6328 OF 15495 ***
    // Wavefunction(s) for diagram number 6328
    // (none)
    // Amplitude(s) for diagram number 6328
    FFV1_0( w_fp[77], w_fp[596], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6329 OF 15495 ***
    // Wavefunction(s) for diagram number 6329
    // (none)
    // Amplitude(s) for diagram number 6329
    FFV1_0( w_fp[77], w_fp[590], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6330 OF 15495 ***
    // Wavefunction(s) for diagram number 6330
    // (none)
    // Amplitude(s) for diagram number 6330
    VVVV1_0( w_fp[523], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[523], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    VVVV4_0( w_fp[523], w_fp[534], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6331 OF 15495 ***
    // Wavefunction(s) for diagram number 6331
    // (none)
    // Amplitude(s) for diagram number 6331
    VVV1_0( w_fp[534], w_fp[7], w_fp[580], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6332 OF 15495 ***
    // Wavefunction(s) for diagram number 6332
    // (none)
    // Amplitude(s) for diagram number 6332
    VVV1_0( w_fp[534], w_fp[5], w_fp[515], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6333 OF 15495 ***
    // Wavefunction(s) for diagram number 6333
    FFV1_1( w_fp[2], w_fp[523], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[595] );
    // Amplitude(s) for diagram number 6333
    FFV1_0( w_fp[29], w_fp[595], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];

    // *** DIAGRAM 6334 OF 15495 ***
    // Wavefunction(s) for diagram number 6334
    // (none)
    // Amplitude(s) for diagram number 6334
    FFV1_0( w_fp[29], w_fp[2], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6335 OF 15495 ***
    // Wavefunction(s) for diagram number 6335
    // (none)
    // Amplitude(s) for diagram number 6335
    FFV1_0( w_fp[77], w_fp[595], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];

    // *** DIAGRAM 6336 OF 15495 ***
    // Wavefunction(s) for diagram number 6336
    // (none)
    // Amplitude(s) for diagram number 6336
    FFV1_0( w_fp[77], w_fp[2], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6337 OF 15495 ***
    // Wavefunction(s) for diagram number 6337
    // (none)
    // Amplitude(s) for diagram number 6337
    VVVV1_0( w_fp[474], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[534], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6338 OF 15495 ***
    // Wavefunction(s) for diagram number 6338
    // (none)
    // Amplitude(s) for diagram number 6338
    VVV1_0( w_fp[534], w_fp[7], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6339 OF 15495 ***
    // Wavefunction(s) for diagram number 6339
    // (none)
    // Amplitude(s) for diagram number 6339
    VVV1_0( w_fp[534], w_fp[4], w_fp[536], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6340 OF 15495 ***
    // Wavefunction(s) for diagram number 6340
    FFV1_1( w_fp[2], w_fp[474], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[594] );
    // Amplitude(s) for diagram number 6340
    FFV1_0( w_fp[94], w_fp[594], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];

    // *** DIAGRAM 6341 OF 15495 ***
    // Wavefunction(s) for diagram number 6341
    // (none)
    // Amplitude(s) for diagram number 6341
    FFV1_0( w_fp[94], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6342 OF 15495 ***
    // Wavefunction(s) for diagram number 6342
    // (none)
    // Amplitude(s) for diagram number 6342
    FFV1_0( w_fp[77], w_fp[594], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];

    // *** DIAGRAM 6343 OF 15495 ***
    // Wavefunction(s) for diagram number 6343
    // (none)
    // Amplitude(s) for diagram number 6343
    FFV1_0( w_fp[77], w_fp[2], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6344 OF 15495 ***
    // Wavefunction(s) for diagram number 6344
    // (none)
    // Amplitude(s) for diagram number 6344
    VVVV1_0( w_fp[547], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVVV3_0( w_fp[547], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    VVVV4_0( w_fp[547], w_fp[534], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6345 OF 15495 ***
    // Wavefunction(s) for diagram number 6345
    // (none)
    // Amplitude(s) for diagram number 6345
    VVV1_0( w_fp[534], w_fp[5], w_fp[522], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6346 OF 15495 ***
    // Wavefunction(s) for diagram number 6346
    // (none)
    // Amplitude(s) for diagram number 6346
    VVV1_0( w_fp[534], w_fp[4], w_fp[519], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6347 OF 15495 ***
    // Wavefunction(s) for diagram number 6347
    FFV1_1( w_fp[2], w_fp[547], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[571] );
    // Amplitude(s) for diagram number 6347
    FFV1_0( w_fp[94], w_fp[571], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 6348 OF 15495 ***
    // Wavefunction(s) for diagram number 6348
    // (none)
    // Amplitude(s) for diagram number 6348
    FFV1_0( w_fp[94], w_fp[2], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6349 OF 15495 ***
    // Wavefunction(s) for diagram number 6349
    // (none)
    // Amplitude(s) for diagram number 6349
    FFV1_0( w_fp[29], w_fp[571], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 6350 OF 15495 ***
    // Wavefunction(s) for diagram number 6350
    // (none)
    // Amplitude(s) for diagram number 6350
    FFV1_0( w_fp[29], w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6351 OF 15495 ***
    // Wavefunction(s) for diagram number 6351
    // (none)
    // Amplitude(s) for diagram number 6351
    VVV1_0( w_fp[588], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0( w_fp[586], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0( w_fp[451], w_fp[534], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6352 OF 15495 ***
    // Wavefunction(s) for diagram number 6352
    // (none)
    // Amplitude(s) for diagram number 6352
    FFV1_0( w_fp[77], w_fp[2], w_fp[588], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[489] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[586], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[495] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[77], w_fp[2], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6353 OF 15495 ***
    // Wavefunction(s) for diagram number 6353
    // (none)
    // Amplitude(s) for diagram number 6353
    VVV1_0( w_fp[520], w_fp[534], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    VVV1_0( w_fp[486], w_fp[534], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    VVV1_0( w_fp[576], w_fp[534], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6354 OF 15495 ***
    // Wavefunction(s) for diagram number 6354
    // (none)
    // Amplitude(s) for diagram number 6354
    FFV1_0( w_fp[29], w_fp[2], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[491] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[29], w_fp[2], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[29], w_fp[2], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6355 OF 15495 ***
    // Wavefunction(s) for diagram number 6355
    // (none)
    // Amplitude(s) for diagram number 6355
    VVV1_0( w_fp[454], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    VVV1_0( w_fp[575], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0( w_fp[589], w_fp[534], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[349] += amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];

    // *** DIAGRAM 6356 OF 15495 ***
    // Wavefunction(s) for diagram number 6356
    // (none)
    // Amplitude(s) for diagram number 6356
    FFV1_0( w_fp[94], w_fp[2], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[94], w_fp[2], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6357 OF 15495 ***
    // Wavefunction(s) for diagram number 6357
    FFV1_2( w_fp[157], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    // Amplitude(s) for diagram number 6357
    FFV1_0( w_fp[495], w_fp[158], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6358 OF 15495 ***
    // Wavefunction(s) for diagram number 6358
    // (none)
    // Amplitude(s) for diagram number 6358
    FFV1_0( w_fp[495], w_fp[163], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6359 OF 15495 ***
    // Wavefunction(s) for diagram number 6359
    FFV1_1( w_fp[156], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[553] );
    // Amplitude(s) for diagram number 6359
    FFV1_0( w_fp[29], w_fp[553], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6360 OF 15495 ***
    // Wavefunction(s) for diagram number 6360
    // (none)
    // Amplitude(s) for diagram number 6360
    FFV1_0( w_fp[77], w_fp[553], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6361 OF 15495 ***
    // Wavefunction(s) for diagram number 6361
    // (none)
    // Amplitude(s) for diagram number 6361
    VVV1_0( w_fp[474], w_fp[497], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6362 OF 15495 ***
    // Wavefunction(s) for diagram number 6362
    // (none)
    // Amplitude(s) for diagram number 6362
    FFV1_0( w_fp[77], w_fp[156], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];

    // *** DIAGRAM 6363 OF 15495 ***
    // Wavefunction(s) for diagram number 6363
    // (none)
    // Amplitude(s) for diagram number 6363
    FFV1_0( w_fp[157], w_fp[163], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[341] += amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];

    // *** DIAGRAM 6364 OF 15495 ***
    // Wavefunction(s) for diagram number 6364
    // (none)
    // Amplitude(s) for diagram number 6364
    VVV1_0( w_fp[547], w_fp[497], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6365 OF 15495 ***
    // Wavefunction(s) for diagram number 6365
    // (none)
    // Amplitude(s) for diagram number 6365
    FFV1_0( w_fp[29], w_fp[156], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];

    // *** DIAGRAM 6366 OF 15495 ***
    // Wavefunction(s) for diagram number 6366
    // (none)
    // Amplitude(s) for diagram number 6366
    FFV1_0( w_fp[157], w_fp[158], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[307] -= amp_sv[0];
    jamp_sv[310] += amp_sv[0];

    // *** DIAGRAM 6367 OF 15495 ***
    // Wavefunction(s) for diagram number 6367
    // (none)
    // Amplitude(s) for diagram number 6367
    FFV1_0( w_fp[29], w_fp[163], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[341] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6368 OF 15495 ***
    // Wavefunction(s) for diagram number 6368
    // (none)
    // Amplitude(s) for diagram number 6368
    FFV1_0( w_fp[77], w_fp[158], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[291] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6369 OF 15495 ***
    // Wavefunction(s) for diagram number 6369
    // (none)
    // Amplitude(s) for diagram number 6369
    FFV1_0( w_fp[157], w_fp[156], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[156], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6370 OF 15495 ***
    // Wavefunction(s) for diagram number 6370
    // (none)
    // Amplitude(s) for diagram number 6370
    FFV1_0( w_fp[495], w_fp[156], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[307] += amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[349] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];

    // *** DIAGRAM 6371 OF 15495 ***
    // Wavefunction(s) for diagram number 6371
    // (none)
    // Amplitude(s) for diagram number 6371
    FFV1_0( w_fp[157], w_fp[553], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];

    // *** DIAGRAM 6372 OF 15495 ***
    // Wavefunction(s) for diagram number 6372
    // (none)
    // Amplitude(s) for diagram number 6372
    FFV1_0( w_fp[157], w_fp[156], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6373 OF 15495 ***
    // Wavefunction(s) for diagram number 6373
    // (none)
    // Amplitude(s) for diagram number 6373
    FFV1_0( w_fp[495], w_fp[190], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6374 OF 15495 ***
    // Wavefunction(s) for diagram number 6374
    // (none)
    // Amplitude(s) for diagram number 6374
    FFV1_0( w_fp[495], w_fp[193], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6375 OF 15495 ***
    // Wavefunction(s) for diagram number 6375
    FFV1_1( w_fp[169], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[535] );
    // Amplitude(s) for diagram number 6375
    FFV1_0( w_fp[94], w_fp[535], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6376 OF 15495 ***
    // Wavefunction(s) for diagram number 6376
    // (none)
    // Amplitude(s) for diagram number 6376
    FFV1_0( w_fp[77], w_fp[535], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6377 OF 15495 ***
    // Wavefunction(s) for diagram number 6377
    // (none)
    // Amplitude(s) for diagram number 6377
    VVV1_0( w_fp[523], w_fp[540], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6378 OF 15495 ***
    // Wavefunction(s) for diagram number 6378
    // (none)
    // Amplitude(s) for diagram number 6378
    FFV1_0( w_fp[77], w_fp[169], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];

    // *** DIAGRAM 6379 OF 15495 ***
    // Wavefunction(s) for diagram number 6379
    // (none)
    // Amplitude(s) for diagram number 6379
    FFV1_0( w_fp[157], w_fp[193], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 6380 OF 15495 ***
    // Wavefunction(s) for diagram number 6380
    // (none)
    // Amplitude(s) for diagram number 6380
    VVV1_0( w_fp[547], w_fp[540], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6381 OF 15495 ***
    // Wavefunction(s) for diagram number 6381
    // (none)
    // Amplitude(s) for diagram number 6381
    FFV1_0( w_fp[94], w_fp[169], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];

    // *** DIAGRAM 6382 OF 15495 ***
    // Wavefunction(s) for diagram number 6382
    // (none)
    // Amplitude(s) for diagram number 6382
    FFV1_0( w_fp[157], w_fp[190], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];

    // *** DIAGRAM 6383 OF 15495 ***
    // Wavefunction(s) for diagram number 6383
    // (none)
    // Amplitude(s) for diagram number 6383
    FFV1_0( w_fp[94], w_fp[193], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6384 OF 15495 ***
    // Wavefunction(s) for diagram number 6384
    // (none)
    // Amplitude(s) for diagram number 6384
    FFV1_0( w_fp[77], w_fp[190], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6385 OF 15495 ***
    // Wavefunction(s) for diagram number 6385
    // (none)
    // Amplitude(s) for diagram number 6385
    FFV1_0( w_fp[157], w_fp[169], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[475] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[169], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6386 OF 15495 ***
    // Wavefunction(s) for diagram number 6386
    // (none)
    // Amplitude(s) for diagram number 6386
    FFV1_0( w_fp[495], w_fp[169], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[427] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];

    // *** DIAGRAM 6387 OF 15495 ***
    // Wavefunction(s) for diagram number 6387
    // (none)
    // Amplitude(s) for diagram number 6387
    FFV1_0( w_fp[157], w_fp[535], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];

    // *** DIAGRAM 6388 OF 15495 ***
    // Wavefunction(s) for diagram number 6388
    // (none)
    // Amplitude(s) for diagram number 6388
    FFV1_0( w_fp[157], w_fp[169], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6389 OF 15495 ***
    // Wavefunction(s) for diagram number 6389
    // (none)
    // Amplitude(s) for diagram number 6389
    FFV1_0( w_fp[495], w_fp[225], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6390 OF 15495 ***
    // Wavefunction(s) for diagram number 6390
    // (none)
    // Amplitude(s) for diagram number 6390
    FFV1_0( w_fp[495], w_fp[226], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6391 OF 15495 ***
    // Wavefunction(s) for diagram number 6391
    FFV1_1( w_fp[215], w_fp[530], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[435] );
    // Amplitude(s) for diagram number 6391
    FFV1_0( w_fp[94], w_fp[435], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6392 OF 15495 ***
    // Wavefunction(s) for diagram number 6392
    // (none)
    // Amplitude(s) for diagram number 6392
    FFV1_0( w_fp[29], w_fp[435], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6393 OF 15495 ***
    // Wavefunction(s) for diagram number 6393
    // (none)
    // Amplitude(s) for diagram number 6393
    VVV1_0( w_fp[523], w_fp[544], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 6394 OF 15495 ***
    // Wavefunction(s) for diagram number 6394
    // (none)
    // Amplitude(s) for diagram number 6394
    FFV1_0( w_fp[29], w_fp[215], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];

    // *** DIAGRAM 6395 OF 15495 ***
    // Wavefunction(s) for diagram number 6395
    // (none)
    // Amplitude(s) for diagram number 6395
    FFV1_0( w_fp[157], w_fp[226], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];

    // *** DIAGRAM 6396 OF 15495 ***
    // Wavefunction(s) for diagram number 6396
    // (none)
    // Amplitude(s) for diagram number 6396
    VVV1_0( w_fp[474], w_fp[544], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6397 OF 15495 ***
    // Wavefunction(s) for diagram number 6397
    // (none)
    // Amplitude(s) for diagram number 6397
    FFV1_0( w_fp[94], w_fp[215], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];

    // *** DIAGRAM 6398 OF 15495 ***
    // Wavefunction(s) for diagram number 6398
    // (none)
    // Amplitude(s) for diagram number 6398
    FFV1_0( w_fp[157], w_fp[225], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];

    // *** DIAGRAM 6399 OF 15495 ***
    // Wavefunction(s) for diagram number 6399
    // (none)
    // Amplitude(s) for diagram number 6399
    FFV1_0( w_fp[94], w_fp[226], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6400 OF 15495 ***
    // Wavefunction(s) for diagram number 6400
    // (none)
    // Amplitude(s) for diagram number 6400
    FFV1_0( w_fp[29], w_fp[225], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 435 );
    storeWf( wfs, w_cx, nevt, 437 );
    storeWf( wfs, w_cx, nevt, 438 );
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 449 );
    storeWf( wfs, w_cx, nevt, 451 );
    storeWf( wfs, w_cx, nevt, 453 );
    storeWf( wfs, w_cx, nevt, 454 );
    storeWf( wfs, w_cx, nevt, 474 );
    storeWf( wfs, w_cx, nevt, 479 );
    storeWf( wfs, w_cx, nevt, 480 );
    storeWf( wfs, w_cx, nevt, 481 );
    storeWf( wfs, w_cx, nevt, 486 );
    storeWf( wfs, w_cx, nevt, 488 );
    storeWf( wfs, w_cx, nevt, 495 );
    storeWf( wfs, w_cx, nevt, 510 );
    storeWf( wfs, w_cx, nevt, 514 );
    storeWf( wfs, w_cx, nevt, 515 );
    storeWf( wfs, w_cx, nevt, 518 );
    storeWf( wfs, w_cx, nevt, 519 );
    storeWf( wfs, w_cx, nevt, 520 );
    storeWf( wfs, w_cx, nevt, 521 );
    storeWf( wfs, w_cx, nevt, 522 );
    storeWf( wfs, w_cx, nevt, 523 );
    storeWf( wfs, w_cx, nevt, 524 );
    storeWf( wfs, w_cx, nevt, 529 );
    storeWf( wfs, w_cx, nevt, 530 );
    storeWf( wfs, w_cx, nevt, 532 );
    storeWf( wfs, w_cx, nevt, 535 );
    storeWf( wfs, w_cx, nevt, 536 );
    storeWf( wfs, w_cx, nevt, 541 );
    storeWf( wfs, w_cx, nevt, 546 );
    storeWf( wfs, w_cx, nevt, 547 );
    storeWf( wfs, w_cx, nevt, 550 );
    storeWf( wfs, w_cx, nevt, 553 );
    storeWf( wfs, w_cx, nevt, 555 );
    storeWf( wfs, w_cx, nevt, 556 );
    storeWf( wfs, w_cx, nevt, 557 );
    storeWf( wfs, w_cx, nevt, 558 );
    storeWf( wfs, w_cx, nevt, 559 );
    storeWf( wfs, w_cx, nevt, 560 );
    storeWf( wfs, w_cx, nevt, 562 );
    storeWf( wfs, w_cx, nevt, 571 );
    storeWf( wfs, w_cx, nevt, 575 );
    storeWf( wfs, w_cx, nevt, 576 );
    storeWf( wfs, w_cx, nevt, 577 );
    storeWf( wfs, w_cx, nevt, 580 );
    storeWf( wfs, w_cx, nevt, 585 );
    storeWf( wfs, w_cx, nevt, 586 );
    storeWf( wfs, w_cx, nevt, 588 );
    storeWf( wfs, w_cx, nevt, 589 );
    storeWf( wfs, w_cx, nevt, 590 );
    storeWf( wfs, w_cx, nevt, 594 );
    storeWf( wfs, w_cx, nevt, 595 );
    storeWf( wfs, w_cx, nevt, 596 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
