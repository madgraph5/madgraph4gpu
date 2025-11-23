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
  diagramgroup13( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 48 );
    retrieveWf( wfs, w_cx, nevt, 49 );
    retrieveWf( wfs, w_cx, nevt, 50 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 97 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 105 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 111 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 119 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 165 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 192 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 240 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 280 );
    retrieveWf( wfs, w_cx, nevt, 282 );
    retrieveWf( wfs, w_cx, nevt, 284 );
    retrieveWf( wfs, w_cx, nevt, 305 );
    retrieveWf( wfs, w_cx, nevt, 345 );
    retrieveWf( wfs, w_cx, nevt, 390 );
    retrieveWf( wfs, w_cx, nevt, 391 );
    retrieveWf( wfs, w_cx, nevt, 400 );
    retrieveWf( wfs, w_cx, nevt, 401 );
    retrieveWf( wfs, w_cx, nevt, 402 );
#endif
#endif

    // *** DIAGRAM 2401 OF 15495 ***
    // Wavefunction(s) for diagram number 2401
    // (none)
    // Amplitude(s) for diagram number 2401
    FFV1_0( w_fp[174], w_fp[95], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[95], w_fp[391], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[95], w_fp[390], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];

    // *** DIAGRAM 2402 OF 15495 ***
    // Wavefunction(s) for diagram number 2402
    // (none)
    // Amplitude(s) for diagram number 2402
    FFV1_0( w_fp[3], w_fp[400], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[401], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[402], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2403 OF 15495 ***
    // Wavefunction(s) for diagram number 2403
    // (none)
    // Amplitude(s) for diagram number 2403
    VVV1_0( w_fp[97], w_fp[113], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    VVV1_0( w_fp[391], w_fp[113], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    VVV1_0( w_fp[390], w_fp[113], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[62] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];

    // *** DIAGRAM 2404 OF 15495 ***
    // Wavefunction(s) for diagram number 2404
    // (none)
    // Amplitude(s) for diagram number 2404
    FFV1_0( w_fp[3], w_fp[127], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[127], w_fp[391], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[127], w_fp[390], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2405 OF 15495 ***
    // Wavefunction(s) for diagram number 2405
    VVVV1P0_1( w_fp[1], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[127] );
    VVVV3P0_1( w_fp[1], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[402] );
    VVVV4P0_1( w_fp[1], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[401] );
    VVV1P0_1( w_fp[127], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[400] );
    VVV1P0_1( w_fp[402], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[403] );
    VVV1P0_1( w_fp[401], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[404] );
    // Amplitude(s) for diagram number 2405
    VVV1_0( w_fp[259], w_fp[400], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[403], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[404], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2406 OF 15495 ***
    // Wavefunction(s) for diagram number 2406
    VVV1P0_1( w_fp[127], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[405] );
    VVV1P0_1( w_fp[402], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[406] );
    VVV1P0_1( w_fp[401], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[407] );
    // Amplitude(s) for diagram number 2406
    VVV1_0( w_fp[259], w_fp[405], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[406], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[407], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2407 OF 15495 ***
    // Wavefunction(s) for diagram number 2407
    // (none)
    // Amplitude(s) for diagram number 2407
    VVVV1_0( w_fp[127], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV3_0( w_fp[127], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV4_0( w_fp[127], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    VVVV1_0( w_fp[402], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVVV3_0( w_fp[402], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVVV4_0( w_fp[402], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    VVVV1_0( w_fp[401], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[401], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[401], w_fp[4], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2408 OF 15495 ***
    // Wavefunction(s) for diagram number 2408
    FFV1_2( w_fp[3], w_fp[127], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[408] );
    FFV1_2( w_fp[3], w_fp[402], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[409] );
    FFV1_2( w_fp[3], w_fp[401], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[410] );
    // Amplitude(s) for diagram number 2408
    FFV1_0( w_fp[408], w_fp[45], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    FFV1_0( w_fp[409], w_fp[45], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    FFV1_0( w_fp[410], w_fp[45], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];

    // *** DIAGRAM 2409 OF 15495 ***
    // Wavefunction(s) for diagram number 2409
    // (none)
    // Amplitude(s) for diagram number 2409
    FFV1_0( w_fp[3], w_fp[45], w_fp[405], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[406], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[407], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2410 OF 15495 ***
    // Wavefunction(s) for diagram number 2410
    // (none)
    // Amplitude(s) for diagram number 2410
    FFV1_0( w_fp[408], w_fp[108], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    FFV1_0( w_fp[409], w_fp[108], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    FFV1_0( w_fp[410], w_fp[108], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];

    // *** DIAGRAM 2411 OF 15495 ***
    // Wavefunction(s) for diagram number 2411
    // (none)
    // Amplitude(s) for diagram number 2411
    FFV1_0( w_fp[3], w_fp[108], w_fp[400], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[108], w_fp[403], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[108], w_fp[404], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2412 OF 15495 ***
    // Wavefunction(s) for diagram number 2412
    FFV1_1( w_fp[10], w_fp[127], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[411] );
    FFV1_1( w_fp[10], w_fp[402], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[412] );
    FFV1_1( w_fp[10], w_fp[401], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[413] );
    // Amplitude(s) for diagram number 2412
    FFV1_0( w_fp[196], w_fp[411], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[412], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[413], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];

    // *** DIAGRAM 2413 OF 15495 ***
    // Wavefunction(s) for diagram number 2413
    // (none)
    // Amplitude(s) for diagram number 2413
    VVV1_0( w_fp[127], w_fp[7], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[402], w_fp[7], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[401], w_fp[7], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2414 OF 15495 ***
    // Wavefunction(s) for diagram number 2414
    // (none)
    // Amplitude(s) for diagram number 2414
    FFV1_0( w_fp[196], w_fp[108], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[108], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[108], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 2415 OF 15495 ***
    // Wavefunction(s) for diagram number 2415
    // (none)
    // Amplitude(s) for diagram number 2415
    FFV1_0( w_fp[179], w_fp[411], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[412], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[60] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[413], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[84] -= amp_sv[0];

    // *** DIAGRAM 2416 OF 15495 ***
    // Wavefunction(s) for diagram number 2416
    // (none)
    // Amplitude(s) for diagram number 2416
    VVV1_0( w_fp[127], w_fp[4], w_fp[284], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[402], w_fp[4], w_fp[284], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[401], w_fp[4], w_fp[284], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2417 OF 15495 ***
    // Wavefunction(s) for diagram number 2417
    // (none)
    // Amplitude(s) for diagram number 2417
    FFV1_0( w_fp[179], w_fp[45], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[45], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[45], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];

    // *** DIAGRAM 2418 OF 15495 ***
    // Wavefunction(s) for diagram number 2418
    // (none)
    // Amplitude(s) for diagram number 2418
    FFV1_0( w_fp[3], w_fp[411], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[412], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[60] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[413], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2419 OF 15495 ***
    // Wavefunction(s) for diagram number 2419
    // (none)
    // Amplitude(s) for diagram number 2419
    VVV1_0( w_fp[127], w_fp[102], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    VVV1_0( w_fp[402], w_fp[102], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[60] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    VVV1_0( w_fp[401], w_fp[102], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];

    // *** DIAGRAM 2420 OF 15495 ***
    // Wavefunction(s) for diagram number 2420
    // (none)
    // Amplitude(s) for diagram number 2420
    FFV1_0( w_fp[3], w_fp[119], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[119], w_fp[402], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[119], w_fp[401], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2421 OF 15495 ***
    // Wavefunction(s) for diagram number 2421
    VVVV1P0_1( w_fp[1], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[119] );
    VVVV3P0_1( w_fp[1], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[413] );
    VVVV4P0_1( w_fp[1], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[412] );
    VVV1P0_1( w_fp[119], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[411] );
    VVV1P0_1( w_fp[413], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[284] );
    VVV1P0_1( w_fp[412], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[414] );
    // Amplitude(s) for diagram number 2421
    VVV1_0( w_fp[259], w_fp[411], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[284], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[414], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2422 OF 15495 ***
    // Wavefunction(s) for diagram number 2422
    VVV1P0_1( w_fp[119], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[415] );
    VVV1P0_1( w_fp[413], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[416] );
    VVV1P0_1( w_fp[412], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[417] );
    // Amplitude(s) for diagram number 2422
    VVV1_0( w_fp[259], w_fp[415], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[416], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[417], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2423 OF 15495 ***
    // Wavefunction(s) for diagram number 2423
    // (none)
    // Amplitude(s) for diagram number 2423
    VVVV1_0( w_fp[119], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    VVVV3_0( w_fp[119], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    VVVV4_0( w_fp[119], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    VVVV1_0( w_fp[413], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    VVVV3_0( w_fp[413], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVVV4_0( w_fp[413], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    VVVV1_0( w_fp[412], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[412], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[412], w_fp[4], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2424 OF 15495 ***
    // Wavefunction(s) for diagram number 2424
    FFV1_2( w_fp[3], w_fp[119], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[418] );
    FFV1_2( w_fp[3], w_fp[413], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[419] );
    FFV1_2( w_fp[3], w_fp[412], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[420] );
    // Amplitude(s) for diagram number 2424
    FFV1_0( w_fp[418], w_fp[45], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    FFV1_0( w_fp[419], w_fp[45], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    FFV1_0( w_fp[420], w_fp[45], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];

    // *** DIAGRAM 2425 OF 15495 ***
    // Wavefunction(s) for diagram number 2425
    // (none)
    // Amplitude(s) for diagram number 2425
    FFV1_0( w_fp[3], w_fp[45], w_fp[415], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[416], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[417], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2426 OF 15495 ***
    // Wavefunction(s) for diagram number 2426
    // (none)
    // Amplitude(s) for diagram number 2426
    FFV1_0( w_fp[418], w_fp[177], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] += amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    FFV1_0( w_fp[419], w_fp[177], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    FFV1_0( w_fp[420], w_fp[177], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];

    // *** DIAGRAM 2427 OF 15495 ***
    // Wavefunction(s) for diagram number 2427
    // (none)
    // Amplitude(s) for diagram number 2427
    FFV1_0( w_fp[3], w_fp[177], w_fp[411], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[177], w_fp[284], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[177], w_fp[414], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 2428 OF 15495 ***
    // Wavefunction(s) for diagram number 2428
    FFV1_1( w_fp[10], w_fp[119], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[421] );
    FFV1_1( w_fp[10], w_fp[413], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[422] );
    FFV1_1( w_fp[10], w_fp[412], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[423] );
    // Amplitude(s) for diagram number 2428
    FFV1_0( w_fp[196], w_fp[421], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[422], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[423], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];

    // *** DIAGRAM 2429 OF 15495 ***
    // Wavefunction(s) for diagram number 2429
    // (none)
    // Amplitude(s) for diagram number 2429
    VVV1_0( w_fp[119], w_fp[6], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[413], w_fp[6], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[412], w_fp[6], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2430 OF 15495 ***
    // Wavefunction(s) for diagram number 2430
    // (none)
    // Amplitude(s) for diagram number 2430
    FFV1_0( w_fp[196], w_fp[177], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[177], w_fp[413], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[177], w_fp[412], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];

    // *** DIAGRAM 2431 OF 15495 ***
    // Wavefunction(s) for diagram number 2431
    // (none)
    // Amplitude(s) for diagram number 2431
    FFV1_0( w_fp[174], w_fp[421], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[422], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[423], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];

    // *** DIAGRAM 2432 OF 15495 ***
    // Wavefunction(s) for diagram number 2432
    // (none)
    // Amplitude(s) for diagram number 2432
    VVV1_0( w_fp[119], w_fp[4], w_fp[282], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[413], w_fp[4], w_fp[282], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[412], w_fp[4], w_fp[282], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2433 OF 15495 ***
    // Wavefunction(s) for diagram number 2433
    // (none)
    // Amplitude(s) for diagram number 2433
    FFV1_0( w_fp[174], w_fp[45], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[45], w_fp[413], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[45], w_fp[412], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];

    // *** DIAGRAM 2434 OF 15495 ***
    // Wavefunction(s) for diagram number 2434
    // (none)
    // Amplitude(s) for diagram number 2434
    FFV1_0( w_fp[3], w_fp[421], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[422], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[423], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2435 OF 15495 ***
    // Wavefunction(s) for diagram number 2435
    // (none)
    // Amplitude(s) for diagram number 2435
    VVV1_0( w_fp[119], w_fp[86], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    VVV1_0( w_fp[413], w_fp[86], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    VVV1_0( w_fp[412], w_fp[86], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];

    // *** DIAGRAM 2436 OF 15495 ***
    // Wavefunction(s) for diagram number 2436
    // (none)
    // Amplitude(s) for diagram number 2436
    FFV1_0( w_fp[3], w_fp[111], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[111], w_fp[413], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[111], w_fp[412], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2437 OF 15495 ***
    // Wavefunction(s) for diagram number 2437
    VVVV1P0_1( w_fp[1], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[111] );
    VVVV3P0_1( w_fp[1], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[423] );
    VVVV4P0_1( w_fp[1], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[422] );
    VVV1P0_1( w_fp[111], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[421] );
    VVV1P0_1( w_fp[423], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[282] );
    VVV1P0_1( w_fp[422], w_fp[4], COUPs[0], 1.0, depCoup, 0., 0., w_fp[424] );
    // Amplitude(s) for diagram number 2437
    VVV1_0( w_fp[259], w_fp[421], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[282], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[424], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];

    // *** DIAGRAM 2438 OF 15495 ***
    // Wavefunction(s) for diagram number 2438
    VVV1P0_1( w_fp[111], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[425] );
    VVV1P0_1( w_fp[423], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[426] );
    VVV1P0_1( w_fp[422], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[427] );
    // Amplitude(s) for diagram number 2438
    VVV1_0( w_fp[259], w_fp[425], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[426], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[427], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];

    // *** DIAGRAM 2439 OF 15495 ***
    // Wavefunction(s) for diagram number 2439
    // (none)
    // Amplitude(s) for diagram number 2439
    VVVV1_0( w_fp[111], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    VVVV3_0( w_fp[111], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    VVVV4_0( w_fp[111], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    VVVV1_0( w_fp[423], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    VVVV3_0( w_fp[423], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    VVVV4_0( w_fp[423], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    VVVV1_0( w_fp[422], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    VVVV3_0( w_fp[422], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVVV4_0( w_fp[422], w_fp[4], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];

    // *** DIAGRAM 2440 OF 15495 ***
    // Wavefunction(s) for diagram number 2440
    FFV1_2( w_fp[3], w_fp[111], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[428] );
    FFV1_2( w_fp[3], w_fp[423], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[429] );
    FFV1_2( w_fp[3], w_fp[422], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[430] );
    // Amplitude(s) for diagram number 2440
    FFV1_0( w_fp[428], w_fp[45], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    FFV1_0( w_fp[429], w_fp[45], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    FFV1_0( w_fp[430], w_fp[45], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];

    // *** DIAGRAM 2441 OF 15495 ***
    // Wavefunction(s) for diagram number 2441
    // (none)
    // Amplitude(s) for diagram number 2441
    FFV1_0( w_fp[3], w_fp[45], w_fp[425], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[426], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[427], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2442 OF 15495 ***
    // Wavefunction(s) for diagram number 2442
    // (none)
    // Amplitude(s) for diagram number 2442
    FFV1_0( w_fp[428], w_fp[95], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    FFV1_0( w_fp[429], w_fp[95], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    FFV1_0( w_fp[430], w_fp[95], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];

    // *** DIAGRAM 2443 OF 15495 ***
    // Wavefunction(s) for diagram number 2443
    // (none)
    // Amplitude(s) for diagram number 2443
    FFV1_0( w_fp[3], w_fp[95], w_fp[421], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[95], w_fp[282], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[95], w_fp[424], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2444 OF 15495 ***
    // Wavefunction(s) for diagram number 2444
    FFV1_1( w_fp[10], w_fp[111], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[431] );
    FFV1_1( w_fp[10], w_fp[423], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[432] );
    FFV1_1( w_fp[10], w_fp[422], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[433] );
    // Amplitude(s) for diagram number 2444
    FFV1_0( w_fp[196], w_fp[431], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[432], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[433], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];

    // *** DIAGRAM 2445 OF 15495 ***
    // Wavefunction(s) for diagram number 2445
    // (none)
    // Amplitude(s) for diagram number 2445
    VVV1_0( w_fp[111], w_fp[5], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[423], w_fp[5], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[422], w_fp[5], w_fp[305], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2446 OF 15495 ***
    // Wavefunction(s) for diagram number 2446
    // (none)
    // Amplitude(s) for diagram number 2446
    FFV1_0( w_fp[196], w_fp[95], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[95], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[95], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[61] += amp_sv[0];
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];

    // *** DIAGRAM 2447 OF 15495 ***
    // Wavefunction(s) for diagram number 2447
    // (none)
    // Amplitude(s) for diagram number 2447
    FFV1_0( w_fp[168], w_fp[431], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[432], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[433], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];

    // *** DIAGRAM 2448 OF 15495 ***
    // Wavefunction(s) for diagram number 2448
    // (none)
    // Amplitude(s) for diagram number 2448
    VVV1_0( w_fp[111], w_fp[4], w_fp[280], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[423], w_fp[4], w_fp[280], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[422], w_fp[4], w_fp[280], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2449 OF 15495 ***
    // Wavefunction(s) for diagram number 2449
    // (none)
    // Amplitude(s) for diagram number 2449
    FFV1_0( w_fp[168], w_fp[45], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[45], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[45], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];

    // *** DIAGRAM 2450 OF 15495 ***
    // Wavefunction(s) for diagram number 2450
    // (none)
    // Amplitude(s) for diagram number 2450
    FFV1_0( w_fp[3], w_fp[431], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[432], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[90] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[433], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2451 OF 15495 ***
    // Wavefunction(s) for diagram number 2451
    // (none)
    // Amplitude(s) for diagram number 2451
    VVV1_0( w_fp[111], w_fp[66], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    VVV1_0( w_fp[423], w_fp[66], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    VVV1_0( w_fp[422], w_fp[66], w_fp[259], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[76] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];

    // *** DIAGRAM 2452 OF 15495 ***
    // Wavefunction(s) for diagram number 2452
    // (none)
    // Amplitude(s) for diagram number 2452
    FFV1_0( w_fp[3], w_fp[91], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[91], w_fp[423], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[91], w_fp[422], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2453 OF 15495 ***
    // Wavefunction(s) for diagram number 2453
    // (none)
    // Amplitude(s) for diagram number 2453
    FFV1_0( w_fp[50], w_fp[345], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    FFV1_0( w_fp[49], w_fp[345], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    FFV1_0( w_fp[48], w_fp[345], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];

    // *** DIAGRAM 2454 OF 15495 ***
    // Wavefunction(s) for diagram number 2454
    // (none)
    // Amplitude(s) for diagram number 2454
    FFV1_0( w_fp[3], w_fp[345], w_fp[165], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[240], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2455 OF 15495 ***
    // Wavefunction(s) for diagram number 2455
    VVV1P0_1( w_fp[1], w_fp[133], COUPs[0], 1.0, depCoup, 0., 0., w_fp[91] );
    VVV1P0_1( w_fp[1], w_fp[134], COUPs[0], 1.0, depCoup, 0., 0., w_fp[433] );
    VVV1P0_1( w_fp[1], w_fp[135], COUPs[0], 1.0, depCoup, 0., 0., w_fp[432] );
    // Amplitude(s) for diagram number 2455
    VVV1_0( w_fp[259], w_fp[91], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[433], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[432], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 2456 OF 15495 ***
    // Wavefunction(s) for diagram number 2456
    // (none)
    // Amplitude(s) for diagram number 2456
    VVV1_0( w_fp[259], w_fp[1], w_fp[165], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[240], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[30], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 2457 OF 15495 ***
    // Wavefunction(s) for diagram number 2457
    // (none)
    // Amplitude(s) for diagram number 2457
    VVVV1_0( w_fp[1], w_fp[133], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[133], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[133], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[134], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[134], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[134], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[135], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[135], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[135], w_fp[7], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 2458 OF 15495 ***
    // Wavefunction(s) for diagram number 2458
    // (none)
    // Amplitude(s) for diagram number 2458
    FFV1_0( w_fp[3], w_fp[108], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[108], w_fp[433], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[113] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[108], w_fp[432], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[105] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[111] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[119] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2459 OF 15495 ***
    // Wavefunction(s) for diagram number 2459
    // (none)
    // Amplitude(s) for diagram number 2459
    FFV1_0( w_fp[50], w_fp[108], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    FFV1_0( w_fp[49], w_fp[108], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    FFV1_0( w_fp[48], w_fp[108], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];

    // *** DIAGRAM 2460 OF 15495 ***
    // Wavefunction(s) for diagram number 2460
    // (none)
    // Amplitude(s) for diagram number 2460
    FFV1_0( w_fp[44], w_fp[345], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    FFV1_0( w_fp[192], w_fp[345], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    FFV1_0( w_fp[75], w_fp[345], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];

    // *** DIAGRAM 2461 OF 15495 ***
    // Wavefunction(s) for diagram number 2461
    // (none)
    // Amplitude(s) for diagram number 2461
    FFV1_0( w_fp[3], w_fp[345], w_fp[162], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2462 OF 15495 ***
    // Wavefunction(s) for diagram number 2462
    VVV1P0_1( w_fp[1], w_fp[139], COUPs[0], 1.0, depCoup, 0., 0., w_fp[108] );
    VVV1P0_1( w_fp[1], w_fp[140], COUPs[0], 1.0, depCoup, 0., 0., w_fp[431] );
    VVV1P0_1( w_fp[1], w_fp[141], COUPs[0], 1.0, depCoup, 0., 0., w_fp[280] );
    // Amplitude(s) for diagram number 2462
    VVV1_0( w_fp[259], w_fp[108], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[431], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[280], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];

    // *** DIAGRAM 2463 OF 15495 ***
    // Wavefunction(s) for diagram number 2463
    // (none)
    // Amplitude(s) for diagram number 2463
    VVV1_0( w_fp[259], w_fp[1], w_fp[162], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[53], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[52], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[105] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];

    // *** DIAGRAM 2464 OF 15495 ***
    // Wavefunction(s) for diagram number 2464
    // (none)
    // Amplitude(s) for diagram number 2464
    VVVV1_0( w_fp[1], w_fp[139], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[139], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[139], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[111] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[140], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[140], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[140], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[141], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[87] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[141], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[141], w_fp[6], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[105] -= amp_sv[0];
    jamp_sv[111] += amp_sv[0];

    // *** DIAGRAM 2465 OF 15495 ***
    // Wavefunction(s) for diagram number 2465
    // (none)
    // Amplitude(s) for diagram number 2465
    FFV1_0( w_fp[3], w_fp[177], w_fp[108], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[177], w_fp[431], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[89] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[177], w_fp[280], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[87] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[95] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2466 OF 15495 ***
    // Wavefunction(s) for diagram number 2466
    // (none)
    // Amplitude(s) for diagram number 2466
    FFV1_0( w_fp[44], w_fp[177], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] += amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[77] += amp_sv[0];
    FFV1_0( w_fp[192], w_fp[177], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    FFV1_0( w_fp[75], w_fp[177], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];

    // *** DIAGRAM 2467 OF 15495 ***
    // Wavefunction(s) for diagram number 2467
    // (none)
    // Amplitude(s) for diagram number 2467
    FFV1_0( w_fp[37], w_fp[345], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    FFV1_0( w_fp[40], w_fp[345], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    FFV1_0( w_fp[12], w_fp[345], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];

    // *** DIAGRAM 2468 OF 15495 ***
    // Wavefunction(s) for diagram number 2468
    // (none)
    // Amplitude(s) for diagram number 2468
    FFV1_0( w_fp[3], w_fp[345], w_fp[159], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[167], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2469 OF 15495 ***
    // Wavefunction(s) for diagram number 2469
    VVV1P0_1( w_fp[1], w_fp[145], COUPs[0], 1.0, depCoup, 0., 0., w_fp[177] );
    VVV1P0_1( w_fp[1], w_fp[146], COUPs[0], 1.0, depCoup, 0., 0., w_fp[305] );
    VVV1P0_1( w_fp[1], w_fp[147], COUPs[0], 1.0, depCoup, 0., 0., w_fp[10] );
    // Amplitude(s) for diagram number 2469
    VVV1_0( w_fp[259], w_fp[177], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[305], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[10], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];

    // *** DIAGRAM 2470 OF 15495 ***
    // Wavefunction(s) for diagram number 2470
    // (none)
    // Amplitude(s) for diagram number 2470
    VVV1_0( w_fp[259], w_fp[1], w_fp[159], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[54], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[167], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];

    // *** DIAGRAM 2471 OF 15495 ***
    // Wavefunction(s) for diagram number 2471
    // (none)
    // Amplitude(s) for diagram number 2471
    VVVV1_0( w_fp[1], w_fp[145], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[145], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[1], w_fp[145], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[146], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[59] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[146], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[1], w_fp[146], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[59] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[147], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[63] -= amp_sv[0];
    jamp_sv[69] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    VVVV3_0( w_fp[1], w_fp[147], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV4_0( w_fp[1], w_fp[147], w_fp[5], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[63] += amp_sv[0];
    jamp_sv[69] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];

    // *** DIAGRAM 2472 OF 15495 ***
    // Wavefunction(s) for diagram number 2472
    // (none)
    // Amplitude(s) for diagram number 2472
    FFV1_0( w_fp[3], w_fp[95], w_fp[177], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[95], w_fp[305], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[59] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[65] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[95], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[63] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[69] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[71] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2473 OF 15495 ***
    // Wavefunction(s) for diagram number 2473
    // (none)
    // Amplitude(s) for diagram number 2473
    FFV1_0( w_fp[37], w_fp[95], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    FFV1_0( w_fp[40], w_fp[95], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    FFV1_0( w_fp[12], w_fp[95], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];

    // *** DIAGRAM 2474 OF 15495 ***
    // Wavefunction(s) for diagram number 2474
    // (none)
    // Amplitude(s) for diagram number 2474
    FFV1_0( w_fp[26], w_fp[345], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    FFV1_0( w_fp[160], w_fp[345], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    FFV1_0( w_fp[105], w_fp[345], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];

    // *** DIAGRAM 2475 OF 15495 ***
    // Wavefunction(s) for diagram number 2475
    // (none)
    // Amplitude(s) for diagram number 2475
    FFV1_0( w_fp[3], w_fp[345], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[345], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2476 OF 15495 ***
    // Wavefunction(s) for diagram number 2476
    VVV1P0_1( w_fp[1], w_fp[151], COUPs[0], 1.0, depCoup, 0., 0., w_fp[345] );
    VVV1P0_1( w_fp[1], w_fp[152], COUPs[0], 1.0, depCoup, 0., 0., w_fp[95] );
    VVV1P0_1( w_fp[1], w_fp[153], COUPs[0], 1.0, depCoup, 0., 0., w_fp[434] );
    // Amplitude(s) for diagram number 2476
    VVV1_0( w_fp[259], w_fp[345], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[95], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    VVV1_0( w_fp[259], w_fp[434], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];

    // *** DIAGRAM 2477 OF 15495 ***
    // Wavefunction(s) for diagram number 2477
    // (none)
    // Amplitude(s) for diagram number 2477
    VVV1_0( w_fp[259], w_fp[1], w_fp[58], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[57], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[71] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[95] += amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    VVV1_0( w_fp[259], w_fp[1], w_fp[38], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[89] -= amp_sv[0];
    jamp_sv[113] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];

    // *** DIAGRAM 2478 OF 15495 ***
    // Wavefunction(s) for diagram number 2478
    // (none)
    // Amplitude(s) for diagram number 2478
    VVVV1_0( w_fp[1], w_fp[4], w_fp[151], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[1], w_fp[4], w_fp[151], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[4], w_fp[151], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[65] += amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[119] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[4], w_fp[152], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[1], w_fp[4], w_fp[152], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[4], w_fp[152], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[71] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[95] -= amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    VVVV1_0( w_fp[1], w_fp[4], w_fp[153], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[1], w_fp[4], w_fp[153], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[88] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    VVVV4_0( w_fp[1], w_fp[4], w_fp[153], w_fp[259], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[65] -= amp_sv[0];
    jamp_sv[89] += amp_sv[0];
    jamp_sv[113] += amp_sv[0];
    jamp_sv[119] -= amp_sv[0];

    // *** DIAGRAM 2479 OF 15495 ***
    // Wavefunction(s) for diagram number 2479
    // (none)
    // Amplitude(s) for diagram number 2479
    FFV1_0( w_fp[3], w_fp[45], w_fp[345], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[95], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[45], w_fp[434], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2480 OF 15495 ***
    // Wavefunction(s) for diagram number 2480
    // (none)
    // Amplitude(s) for diagram number 2480
    FFV1_0( w_fp[26], w_fp[45], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    FFV1_0( w_fp[160], w_fp[45], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    FFV1_0( w_fp[105], w_fp[45], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];

    // *** DIAGRAM 2481 OF 15495 ***
    // Wavefunction(s) for diagram number 2481
    FFV1_2( w_fp[3], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[45] );
    FFV1_1( w_fp[2], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[259] );
    FFV1_2( w_fp[45], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[435] );
    FFV1_1( w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[436] );
    FFV1_2( w_fp[435], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[437] );
    // Amplitude(s) for diagram number 2481
    FFV1_0( w_fp[437], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] -= amp_sv[0];

    // *** DIAGRAM 2482 OF 15495 ***
    // Wavefunction(s) for diagram number 2482
    FFV1_2( w_fp[435], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[438] );
    // Amplitude(s) for diagram number 2482
    FFV1_0( w_fp[438], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] -= amp_sv[0];

    // *** DIAGRAM 2483 OF 15495 ***
    // Wavefunction(s) for diagram number 2483
    FFV1_1( w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[439] );
    FFV1_2( w_fp[435], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[440] );
    // Amplitude(s) for diagram number 2483
    FFV1_0( w_fp[440], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] -= amp_sv[0];

    // *** DIAGRAM 2484 OF 15495 ***
    // Wavefunction(s) for diagram number 2484
    // (none)
    // Amplitude(s) for diagram number 2484
    FFV1_0( w_fp[438], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] -= amp_sv[0];

    // *** DIAGRAM 2485 OF 15495 ***
    // Wavefunction(s) for diagram number 2485
    FFV1_1( w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[441] );
    // Amplitude(s) for diagram number 2485
    FFV1_0( w_fp[440], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 2486 OF 15495 ***
    // Wavefunction(s) for diagram number 2486
    // (none)
    // Amplitude(s) for diagram number 2486
    FFV1_0( w_fp[437], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] -= amp_sv[0];

    // *** DIAGRAM 2487 OF 15495 ***
    // Wavefunction(s) for diagram number 2487
    FFV1_2( w_fp[45], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[442] );
    FFV1_1( w_fp[259], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[443] );
    FFV1_2( w_fp[442], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[444] );
    // Amplitude(s) for diagram number 2487
    FFV1_0( w_fp[444], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= amp_sv[0];

    // *** DIAGRAM 2488 OF 15495 ***
    // Wavefunction(s) for diagram number 2488
    FFV1_2( w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[445] );
    // Amplitude(s) for diagram number 2488
    FFV1_0( w_fp[445], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];

    // *** DIAGRAM 2489 OF 15495 ***
    // Wavefunction(s) for diagram number 2489
    FFV1_2( w_fp[442], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[446] );
    // Amplitude(s) for diagram number 2489
    FFV1_0( w_fp[446], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] -= amp_sv[0];

    // *** DIAGRAM 2490 OF 15495 ***
    // Wavefunction(s) for diagram number 2490
    // (none)
    // Amplitude(s) for diagram number 2490
    FFV1_0( w_fp[445], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] -= amp_sv[0];

    // *** DIAGRAM 2491 OF 15495 ***
    // Wavefunction(s) for diagram number 2491
    // (none)
    // Amplitude(s) for diagram number 2491
    FFV1_0( w_fp[446], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] -= amp_sv[0];

    // *** DIAGRAM 2492 OF 15495 ***
    // Wavefunction(s) for diagram number 2492
    // (none)
    // Amplitude(s) for diagram number 2492
    FFV1_0( w_fp[444], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] -= amp_sv[0];

    // *** DIAGRAM 2493 OF 15495 ***
    // Wavefunction(s) for diagram number 2493
    FFV1_2( w_fp[45], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[447] );
    FFV1_2( w_fp[447], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[448] );
    // Amplitude(s) for diagram number 2493
    FFV1_0( w_fp[448], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];

    // *** DIAGRAM 2494 OF 15495 ***
    // Wavefunction(s) for diagram number 2494
    FFV1_2( w_fp[447], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[449] );
    // Amplitude(s) for diagram number 2494
    FFV1_0( w_fp[449], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];

    // *** DIAGRAM 2495 OF 15495 ***
    // Wavefunction(s) for diagram number 2495
    FFV1_2( w_fp[447], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[450] );
    // Amplitude(s) for diagram number 2495
    FFV1_0( w_fp[450], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] -= amp_sv[0];

    // *** DIAGRAM 2496 OF 15495 ***
    // Wavefunction(s) for diagram number 2496
    // (none)
    // Amplitude(s) for diagram number 2496
    FFV1_0( w_fp[449], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] -= amp_sv[0];

    // *** DIAGRAM 2497 OF 15495 ***
    // Wavefunction(s) for diagram number 2497
    // (none)
    // Amplitude(s) for diagram number 2497
    FFV1_0( w_fp[450], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] -= amp_sv[0];

    // *** DIAGRAM 2498 OF 15495 ***
    // Wavefunction(s) for diagram number 2498
    // (none)
    // Amplitude(s) for diagram number 2498
    FFV1_0( w_fp[448], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] -= amp_sv[0];

    // *** DIAGRAM 2499 OF 15495 ***
    // Wavefunction(s) for diagram number 2499
    FFV1_2( w_fp[45], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[451] );
    FFV1_2( w_fp[451], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[452] );
    // Amplitude(s) for diagram number 2499
    FFV1_0( w_fp[452], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];

    // *** DIAGRAM 2500 OF 15495 ***
    // Wavefunction(s) for diagram number 2500
    FFV1_2( w_fp[451], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[453] );
    // Amplitude(s) for diagram number 2500
    FFV1_0( w_fp[453], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];

    // *** DIAGRAM 2501 OF 15495 ***
    // Wavefunction(s) for diagram number 2501
    FFV1_2( w_fp[451], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[454] );
    // Amplitude(s) for diagram number 2501
    FFV1_0( w_fp[454], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] -= amp_sv[0];

    // *** DIAGRAM 2502 OF 15495 ***
    // Wavefunction(s) for diagram number 2502
    // (none)
    // Amplitude(s) for diagram number 2502
    FFV1_0( w_fp[453], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= amp_sv[0];

    // *** DIAGRAM 2503 OF 15495 ***
    // Wavefunction(s) for diagram number 2503
    // (none)
    // Amplitude(s) for diagram number 2503
    FFV1_0( w_fp[454], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] -= amp_sv[0];

    // *** DIAGRAM 2504 OF 15495 ***
    // Wavefunction(s) for diagram number 2504
    // (none)
    // Amplitude(s) for diagram number 2504
    FFV1_0( w_fp[452], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] -= amp_sv[0];

    // *** DIAGRAM 2505 OF 15495 ***
    // Wavefunction(s) for diagram number 2505
    FFV1P0_3( w_fp[45], w_fp[259], COUPs[1], 1.0, depCoup, 0., 0., w_fp[455] );
    // Amplitude(s) for diagram number 2505
    VVV1_0( w_fp[455], w_fp[68], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2506 OF 15495 ***
    // Wavefunction(s) for diagram number 2506
    // (none)
    // Amplitude(s) for diagram number 2506
    VVV1_0( w_fp[455], w_fp[69], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2507 OF 15495 ***
    // Wavefunction(s) for diagram number 2507
    // (none)
    // Amplitude(s) for diagram number 2507
    VVVV1_0( w_fp[66], w_fp[6], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[66], w_fp[6], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[66], w_fp[6], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2508 OF 15495 ***
    // Wavefunction(s) for diagram number 2508
    FFV1_2( w_fp[45], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[456] );
    // Amplitude(s) for diagram number 2508
    FFV1_0( w_fp[456], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2509 OF 15495 ***
    // Wavefunction(s) for diagram number 2509
    // (none)
    // Amplitude(s) for diagram number 2509
    FFV1_0( w_fp[456], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2510 OF 15495 ***
    // Wavefunction(s) for diagram number 2510
    FFV1_1( w_fp[259], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[457] );
    // Amplitude(s) for diagram number 2510
    FFV1_0( w_fp[447], w_fp[457], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2511 OF 15495 ***
    // Wavefunction(s) for diagram number 2511
    // (none)
    // Amplitude(s) for diagram number 2511
    FFV1_0( w_fp[447], w_fp[441], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2512 OF 15495 ***
    // Wavefunction(s) for diagram number 2512
    // (none)
    // Amplitude(s) for diagram number 2512
    FFV1_0( w_fp[447], w_fp[259], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 2513 OF 15495 ***
    // Wavefunction(s) for diagram number 2513
    // (none)
    // Amplitude(s) for diagram number 2513
    FFV1_0( w_fp[451], w_fp[457], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2514 OF 15495 ***
    // Wavefunction(s) for diagram number 2514
    // (none)
    // Amplitude(s) for diagram number 2514
    FFV1_0( w_fp[451], w_fp[439], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2515 OF 15495 ***
    // Wavefunction(s) for diagram number 2515
    // (none)
    // Amplitude(s) for diagram number 2515
    FFV1_0( w_fp[451], w_fp[259], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];

    // *** DIAGRAM 2516 OF 15495 ***
    // Wavefunction(s) for diagram number 2516
    // (none)
    // Amplitude(s) for diagram number 2516
    FFV1_0( w_fp[45], w_fp[439], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];

    // *** DIAGRAM 2517 OF 15495 ***
    // Wavefunction(s) for diagram number 2517
    // (none)
    // Amplitude(s) for diagram number 2517
    FFV1_0( w_fp[45], w_fp[441], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 2518 OF 15495 ***
    // Wavefunction(s) for diagram number 2518
    // (none)
    // Amplitude(s) for diagram number 2518
    VVV1_0( w_fp[66], w_fp[84], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2519 OF 15495 ***
    // Wavefunction(s) for diagram number 2519
    // (none)
    // Amplitude(s) for diagram number 2519
    FFV1_0( w_fp[456], w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 2520 OF 15495 ***
    // Wavefunction(s) for diagram number 2520
    FFV1_2( w_fp[45], w_fp[84], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[458] );
    // Amplitude(s) for diagram number 2520
    FFV1_0( w_fp[458], w_fp[259], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];

    // *** DIAGRAM 2521 OF 15495 ***
    // Wavefunction(s) for diagram number 2521
    // (none)
    // Amplitude(s) for diagram number 2521
    VVV1_0( w_fp[455], w_fp[87], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2522 OF 15495 ***
    // Wavefunction(s) for diagram number 2522
    // (none)
    // Amplitude(s) for diagram number 2522
    VVV1_0( w_fp[455], w_fp[88], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2523 OF 15495 ***
    // Wavefunction(s) for diagram number 2523
    // (none)
    // Amplitude(s) for diagram number 2523
    VVVV1_0( w_fp[86], w_fp[5], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[86], w_fp[5], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[86], w_fp[5], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2524 OF 15495 ***
    // Wavefunction(s) for diagram number 2524
    FFV1_2( w_fp[45], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[459] );
    // Amplitude(s) for diagram number 2524
    FFV1_0( w_fp[459], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2525 OF 15495 ***
    // Wavefunction(s) for diagram number 2525
    // (none)
    // Amplitude(s) for diagram number 2525
    FFV1_0( w_fp[459], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2526 OF 15495 ***
    // Wavefunction(s) for diagram number 2526
    FFV1_1( w_fp[259], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[460] );
    // Amplitude(s) for diagram number 2526
    FFV1_0( w_fp[442], w_fp[460], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2527 OF 15495 ***
    // Wavefunction(s) for diagram number 2527
    // (none)
    // Amplitude(s) for diagram number 2527
    FFV1_0( w_fp[442], w_fp[441], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2528 OF 15495 ***
    // Wavefunction(s) for diagram number 2528
    // (none)
    // Amplitude(s) for diagram number 2528
    FFV1_0( w_fp[442], w_fp[259], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 2529 OF 15495 ***
    // Wavefunction(s) for diagram number 2529
    // (none)
    // Amplitude(s) for diagram number 2529
    FFV1_0( w_fp[451], w_fp[460], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2530 OF 15495 ***
    // Wavefunction(s) for diagram number 2530
    // (none)
    // Amplitude(s) for diagram number 2530
    FFV1_0( w_fp[451], w_fp[436], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2531 OF 15495 ***
    // Wavefunction(s) for diagram number 2531
    // (none)
    // Amplitude(s) for diagram number 2531
    FFV1_0( w_fp[451], w_fp[259], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];

    // *** DIAGRAM 2532 OF 15495 ***
    // Wavefunction(s) for diagram number 2532
    // (none)
    // Amplitude(s) for diagram number 2532
    FFV1_0( w_fp[45], w_fp[436], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];

    // *** DIAGRAM 2533 OF 15495 ***
    // Wavefunction(s) for diagram number 2533
    // (none)
    // Amplitude(s) for diagram number 2533
    FFV1_0( w_fp[45], w_fp[441], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];

    // *** DIAGRAM 2534 OF 15495 ***
    // Wavefunction(s) for diagram number 2534
    // (none)
    // Amplitude(s) for diagram number 2534
    VVV1_0( w_fp[86], w_fp[100], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2535 OF 15495 ***
    // Wavefunction(s) for diagram number 2535
    // (none)
    // Amplitude(s) for diagram number 2535
    FFV1_0( w_fp[459], w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];

    // *** DIAGRAM 2536 OF 15495 ***
    // Wavefunction(s) for diagram number 2536
    FFV1_2( w_fp[45], w_fp[100], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[461] );
    // Amplitude(s) for diagram number 2536
    FFV1_0( w_fp[461], w_fp[259], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];

    // *** DIAGRAM 2537 OF 15495 ***
    // Wavefunction(s) for diagram number 2537
    // (none)
    // Amplitude(s) for diagram number 2537
    VVV1_0( w_fp[455], w_fp[103], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2538 OF 15495 ***
    // Wavefunction(s) for diagram number 2538
    // (none)
    // Amplitude(s) for diagram number 2538
    VVV1_0( w_fp[455], w_fp[104], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2539 OF 15495 ***
    // Wavefunction(s) for diagram number 2539
    // (none)
    // Amplitude(s) for diagram number 2539
    VVVV1_0( w_fp[102], w_fp[5], w_fp[6], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[102], w_fp[5], w_fp[6], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[102], w_fp[5], w_fp[6], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2540 OF 15495 ***
    // Wavefunction(s) for diagram number 2540
    FFV1_2( w_fp[45], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[462] );
    // Amplitude(s) for diagram number 2540
    FFV1_0( w_fp[462], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2541 OF 15495 ***
    // Wavefunction(s) for diagram number 2541
    // (none)
    // Amplitude(s) for diagram number 2541
    FFV1_0( w_fp[462], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2542 OF 15495 ***
    // Wavefunction(s) for diagram number 2542
    FFV1_1( w_fp[259], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[463] );
    // Amplitude(s) for diagram number 2542
    FFV1_0( w_fp[442], w_fp[463], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2543 OF 15495 ***
    // Wavefunction(s) for diagram number 2543
    // (none)
    // Amplitude(s) for diagram number 2543
    FFV1_0( w_fp[442], w_fp[439], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2544 OF 15495 ***
    // Wavefunction(s) for diagram number 2544
    // (none)
    // Amplitude(s) for diagram number 2544
    FFV1_0( w_fp[442], w_fp[259], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];

    // *** DIAGRAM 2545 OF 15495 ***
    // Wavefunction(s) for diagram number 2545
    // (none)
    // Amplitude(s) for diagram number 2545
    FFV1_0( w_fp[447], w_fp[463], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2546 OF 15495 ***
    // Wavefunction(s) for diagram number 2546
    // (none)
    // Amplitude(s) for diagram number 2546
    FFV1_0( w_fp[447], w_fp[436], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2547 OF 15495 ***
    // Wavefunction(s) for diagram number 2547
    // (none)
    // Amplitude(s) for diagram number 2547
    FFV1_0( w_fp[447], w_fp[259], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];

    // *** DIAGRAM 2548 OF 15495 ***
    // Wavefunction(s) for diagram number 2548
    // (none)
    // Amplitude(s) for diagram number 2548
    FFV1_0( w_fp[45], w_fp[436], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];

    // *** DIAGRAM 2549 OF 15495 ***
    // Wavefunction(s) for diagram number 2549
    // (none)
    // Amplitude(s) for diagram number 2549
    FFV1_0( w_fp[45], w_fp[439], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];

    // *** DIAGRAM 2550 OF 15495 ***
    // Wavefunction(s) for diagram number 2550
    // (none)
    // Amplitude(s) for diagram number 2550
    VVV1_0( w_fp[102], w_fp[113], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2551 OF 15495 ***
    // Wavefunction(s) for diagram number 2551
    // (none)
    // Amplitude(s) for diagram number 2551
    FFV1_0( w_fp[462], w_fp[259], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];

    // *** DIAGRAM 2552 OF 15495 ***
    // Wavefunction(s) for diagram number 2552
    FFV1_2( w_fp[45], w_fp[113], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[464] );
    // Amplitude(s) for diagram number 2552
    FFV1_0( w_fp[464], w_fp[259], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];

    // *** DIAGRAM 2553 OF 15495 ***
    // Wavefunction(s) for diagram number 2553
    // (none)
    // Amplitude(s) for diagram number 2553
    VVV1_0( w_fp[455], w_fp[115], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2554 OF 15495 ***
    // Wavefunction(s) for diagram number 2554
    // (none)
    // Amplitude(s) for diagram number 2554
    VVV1_0( w_fp[455], w_fp[4], w_fp[116], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2555 OF 15495 ***
    // Wavefunction(s) for diagram number 2555
    // (none)
    // Amplitude(s) for diagram number 2555
    VVVV1_0( w_fp[4], w_fp[113], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[4], w_fp[113], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[4], w_fp[113], w_fp[7], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2556 OF 15495 ***
    // Wavefunction(s) for diagram number 2556
    FFV1_1( w_fp[259], w_fp[113], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[465] );
    // Amplitude(s) for diagram number 2556
    FFV1_0( w_fp[435], w_fp[465], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2557 OF 15495 ***
    // Wavefunction(s) for diagram number 2557
    // (none)
    // Amplitude(s) for diagram number 2557
    FFV1_0( w_fp[435], w_fp[441], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2558 OF 15495 ***
    // Wavefunction(s) for diagram number 2558
    // (none)
    // Amplitude(s) for diagram number 2558
    FFV1_0( w_fp[435], w_fp[259], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 2559 OF 15495 ***
    // Wavefunction(s) for diagram number 2559
    // (none)
    // Amplitude(s) for diagram number 2559
    FFV1_0( w_fp[464], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2560 OF 15495 ***
    // Wavefunction(s) for diagram number 2560
    // (none)
    // Amplitude(s) for diagram number 2560
    FFV1_0( w_fp[464], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2561 OF 15495 ***
    // Wavefunction(s) for diagram number 2561
    // (none)
    // Amplitude(s) for diagram number 2561
    FFV1_0( w_fp[451], w_fp[443], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2562 OF 15495 ***
    // Wavefunction(s) for diagram number 2562
    // (none)
    // Amplitude(s) for diagram number 2562
    FFV1_0( w_fp[451], w_fp[465], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2563 OF 15495 ***
    // Wavefunction(s) for diagram number 2563
    // (none)
    // Amplitude(s) for diagram number 2563
    FFV1_0( w_fp[451], w_fp[259], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];

    // *** DIAGRAM 2564 OF 15495 ***
    // Wavefunction(s) for diagram number 2564
    // (none)
    // Amplitude(s) for diagram number 2564
    FFV1_0( w_fp[45], w_fp[443], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];

    // *** DIAGRAM 2565 OF 15495 ***
    // Wavefunction(s) for diagram number 2565
    // (none)
    // Amplitude(s) for diagram number 2565
    FFV1_0( w_fp[45], w_fp[441], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 2566 OF 15495 ***
    // Wavefunction(s) for diagram number 2566
    // (none)
    // Amplitude(s) for diagram number 2566
    VVV1_0( w_fp[455], w_fp[124], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2567 OF 15495 ***
    // Wavefunction(s) for diagram number 2567
    // (none)
    // Amplitude(s) for diagram number 2567
    VVV1_0( w_fp[455], w_fp[4], w_fp[125], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2568 OF 15495 ***
    // Wavefunction(s) for diagram number 2568
    // (none)
    // Amplitude(s) for diagram number 2568
    VVVV1_0( w_fp[4], w_fp[100], w_fp[6], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[4], w_fp[100], w_fp[6], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[4], w_fp[100], w_fp[6], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2569 OF 15495 ***
    // Wavefunction(s) for diagram number 2569
    FFV1_1( w_fp[259], w_fp[100], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[466] );
    // Amplitude(s) for diagram number 2569
    FFV1_0( w_fp[435], w_fp[466], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2570 OF 15495 ***
    // Wavefunction(s) for diagram number 2570
    // (none)
    // Amplitude(s) for diagram number 2570
    FFV1_0( w_fp[435], w_fp[439], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[209] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2571 OF 15495 ***
    // Wavefunction(s) for diagram number 2571
    // (none)
    // Amplitude(s) for diagram number 2571
    FFV1_0( w_fp[435], w_fp[259], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[191] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];

    // *** DIAGRAM 2572 OF 15495 ***
    // Wavefunction(s) for diagram number 2572
    // (none)
    // Amplitude(s) for diagram number 2572
    FFV1_0( w_fp[461], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2573 OF 15495 ***
    // Wavefunction(s) for diagram number 2573
    // (none)
    // Amplitude(s) for diagram number 2573
    FFV1_0( w_fp[461], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2574 OF 15495 ***
    // Wavefunction(s) for diagram number 2574
    // (none)
    // Amplitude(s) for diagram number 2574
    FFV1_0( w_fp[447], w_fp[443], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2575 OF 15495 ***
    // Wavefunction(s) for diagram number 2575
    // (none)
    // Amplitude(s) for diagram number 2575
    FFV1_0( w_fp[447], w_fp[466], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2576 OF 15495 ***
    // Wavefunction(s) for diagram number 2576
    // (none)
    // Amplitude(s) for diagram number 2576
    FFV1_0( w_fp[447], w_fp[259], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 2577 OF 15495 ***
    // Wavefunction(s) for diagram number 2577
    // (none)
    // Amplitude(s) for diagram number 2577
    FFV1_0( w_fp[45], w_fp[443], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[161] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];

    // *** DIAGRAM 2578 OF 15495 ***
    // Wavefunction(s) for diagram number 2578
    // (none)
    // Amplitude(s) for diagram number 2578
    FFV1_0( w_fp[45], w_fp[439], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];

    // *** DIAGRAM 2579 OF 15495 ***
    // Wavefunction(s) for diagram number 2579
    // (none)
    // Amplitude(s) for diagram number 2579
    VVV1_0( w_fp[455], w_fp[130], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2580 OF 15495 ***
    // Wavefunction(s) for diagram number 2580
    // (none)
    // Amplitude(s) for diagram number 2580
    VVV1_0( w_fp[455], w_fp[4], w_fp[131], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2581 OF 15495 ***
    // Wavefunction(s) for diagram number 2581
    // (none)
    // Amplitude(s) for diagram number 2581
    VVVV1_0( w_fp[4], w_fp[5], w_fp[84], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0( w_fp[4], w_fp[5], w_fp[84], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0( w_fp[4], w_fp[5], w_fp[84], w_fp[455], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2582 OF 15495 ***
    // Wavefunction(s) for diagram number 2582
    // (none)
    // Amplitude(s) for diagram number 2582
    FFV1_0( w_fp[435], w_fp[436], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2583 OF 15495 ***
    // Wavefunction(s) for diagram number 2583
    FFV1_1( w_fp[259], w_fp[84], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[467] );
    // Amplitude(s) for diagram number 2583
    FFV1_0( w_fp[435], w_fp[467], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2584 OF 15495 ***
    // Wavefunction(s) for diagram number 2584
    // (none)
    // Amplitude(s) for diagram number 2584
    FFV1_0( w_fp[435], w_fp[259], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[185] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 2585 OF 15495 ***
    // Wavefunction(s) for diagram number 2585
    // (none)
    // Amplitude(s) for diagram number 2585
    FFV1_0( w_fp[442], w_fp[443], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2586 OF 15495 ***
    // Wavefunction(s) for diagram number 2586
    // (none)
    // Amplitude(s) for diagram number 2586
    FFV1_0( w_fp[442], w_fp[467], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2587 OF 15495 ***
    // Wavefunction(s) for diagram number 2587
    // (none)
    // Amplitude(s) for diagram number 2587
    FFV1_0( w_fp[442], w_fp[259], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];

    // *** DIAGRAM 2588 OF 15495 ***
    // Wavefunction(s) for diagram number 2588
    // (none)
    // Amplitude(s) for diagram number 2588
    FFV1_0( w_fp[458], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2589 OF 15495 ***
    // Wavefunction(s) for diagram number 2589
    // (none)
    // Amplitude(s) for diagram number 2589
    FFV1_0( w_fp[458], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2590 OF 15495 ***
    // Wavefunction(s) for diagram number 2590
    // (none)
    // Amplitude(s) for diagram number 2590
    FFV1_0( w_fp[45], w_fp[443], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[167] += amp_sv[0];

    // *** DIAGRAM 2591 OF 15495 ***
    // Wavefunction(s) for diagram number 2591
    // (none)
    // Amplitude(s) for diagram number 2591
    FFV1_0( w_fp[45], w_fp[436], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];

    // *** DIAGRAM 2592 OF 15495 ***
    // Wavefunction(s) for diagram number 2592
    // (none)
    // Amplitude(s) for diagram number 2592
    VVV1_0( w_fp[133], w_fp[7], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[134], w_fp[7], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[135], w_fp[7], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2593 OF 15495 ***
    // Wavefunction(s) for diagram number 2593
    FFV1_2( w_fp[45], w_fp[133], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[468] );
    FFV1_2( w_fp[45], w_fp[134], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[469] );
    FFV1_2( w_fp[45], w_fp[135], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[470] );
    // Amplitude(s) for diagram number 2593
    FFV1_0( w_fp[468], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] += amp_sv[0];
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    FFV1_0( w_fp[469], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[227] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    FFV1_0( w_fp[470], w_fp[259], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 2594 OF 15495 ***
    // Wavefunction(s) for diagram number 2594
    // (none)
    // Amplitude(s) for diagram number 2594
    FFV1_0( w_fp[451], w_fp[259], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] += amp_sv[0];
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[259], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[159] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    FFV1_0( w_fp[451], w_fp[259], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];

    // *** DIAGRAM 2595 OF 15495 ***
    // Wavefunction(s) for diagram number 2595
    // (none)
    // Amplitude(s) for diagram number 2595
    VVV1_0( w_fp[139], w_fp[6], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[140], w_fp[6], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[141], w_fp[6], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2596 OF 15495 ***
    // Wavefunction(s) for diagram number 2596
    FFV1_2( w_fp[45], w_fp[139], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[471] );
    FFV1_2( w_fp[45], w_fp[140], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[472] );
    FFV1_2( w_fp[45], w_fp[141], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[473] );
    // Amplitude(s) for diagram number 2596
    FFV1_0( w_fp[471], w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] += amp_sv[0];
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    FFV1_0( w_fp[472], w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[203] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    FFV1_0( w_fp[473], w_fp[259], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];

    // *** DIAGRAM 2597 OF 15495 ***
    // Wavefunction(s) for diagram number 2597
    // (none)
    // Amplitude(s) for diagram number 2597
    FFV1_0( w_fp[447], w_fp[259], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] += amp_sv[0];
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[259], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[165] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    FFV1_0( w_fp[447], w_fp[259], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];

    // *** DIAGRAM 2598 OF 15495 ***
    // Wavefunction(s) for diagram number 2598
    // (none)
    // Amplitude(s) for diagram number 2598
    VVV1_0( w_fp[145], w_fp[5], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[146], w_fp[5], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    VVV1_0( w_fp[147], w_fp[5], w_fp[455], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[203] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[227] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 2599 OF 15495 ***
    // Wavefunction(s) for diagram number 2599
    FFV1_2( w_fp[45], w_fp[145], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[474] );
    FFV1_2( w_fp[45], w_fp[146], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[475] );
    FFV1_2( w_fp[45], w_fp[147], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[476] );
    // Amplitude(s) for diagram number 2599
    FFV1_0( w_fp[474], w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    FFV1_0( w_fp[475], w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    FFV1_0( w_fp[476], w_fp[259], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];

    // *** DIAGRAM 2600 OF 15495 ***
    // Wavefunction(s) for diagram number 2600
    // (none)
    // Amplitude(s) for diagram number 2600
    FFV1_0( w_fp[442], w_fp[259], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] += amp_sv[0];
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[259], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[167] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    FFV1_0( w_fp[442], w_fp[259], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[161] -= amp_sv[0];
    jamp_sv[203] += amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 10 );
    storeWf( wfs, w_cx, nevt, 45 );
    storeWf( wfs, w_cx, nevt, 91 );
    storeWf( wfs, w_cx, nevt, 95 );
    storeWf( wfs, w_cx, nevt, 108 );
    storeWf( wfs, w_cx, nevt, 111 );
    storeWf( wfs, w_cx, nevt, 119 );
    storeWf( wfs, w_cx, nevt, 127 );
    storeWf( wfs, w_cx, nevt, 177 );
    storeWf( wfs, w_cx, nevt, 259 );
    storeWf( wfs, w_cx, nevt, 280 );
    storeWf( wfs, w_cx, nevt, 282 );
    storeWf( wfs, w_cx, nevt, 284 );
    storeWf( wfs, w_cx, nevt, 305 );
    storeWf( wfs, w_cx, nevt, 345 );
    storeWf( wfs, w_cx, nevt, 400 );
    storeWf( wfs, w_cx, nevt, 401 );
    storeWf( wfs, w_cx, nevt, 402 );
    storeWf( wfs, w_cx, nevt, 403 );
    storeWf( wfs, w_cx, nevt, 404 );
    storeWf( wfs, w_cx, nevt, 405 );
    storeWf( wfs, w_cx, nevt, 406 );
    storeWf( wfs, w_cx, nevt, 407 );
    storeWf( wfs, w_cx, nevt, 408 );
    storeWf( wfs, w_cx, nevt, 409 );
    storeWf( wfs, w_cx, nevt, 410 );
    storeWf( wfs, w_cx, nevt, 411 );
    storeWf( wfs, w_cx, nevt, 412 );
    storeWf( wfs, w_cx, nevt, 413 );
    storeWf( wfs, w_cx, nevt, 414 );
    storeWf( wfs, w_cx, nevt, 415 );
    storeWf( wfs, w_cx, nevt, 416 );
    storeWf( wfs, w_cx, nevt, 417 );
    storeWf( wfs, w_cx, nevt, 418 );
    storeWf( wfs, w_cx, nevt, 419 );
    storeWf( wfs, w_cx, nevt, 420 );
    storeWf( wfs, w_cx, nevt, 421 );
    storeWf( wfs, w_cx, nevt, 422 );
    storeWf( wfs, w_cx, nevt, 423 );
    storeWf( wfs, w_cx, nevt, 424 );
    storeWf( wfs, w_cx, nevt, 425 );
    storeWf( wfs, w_cx, nevt, 426 );
    storeWf( wfs, w_cx, nevt, 427 );
    storeWf( wfs, w_cx, nevt, 428 );
    storeWf( wfs, w_cx, nevt, 429 );
    storeWf( wfs, w_cx, nevt, 430 );
    storeWf( wfs, w_cx, nevt, 431 );
    storeWf( wfs, w_cx, nevt, 432 );
    storeWf( wfs, w_cx, nevt, 433 );
    storeWf( wfs, w_cx, nevt, 434 );
    storeWf( wfs, w_cx, nevt, 435 );
    storeWf( wfs, w_cx, nevt, 436 );
    storeWf( wfs, w_cx, nevt, 437 );
    storeWf( wfs, w_cx, nevt, 438 );
    storeWf( wfs, w_cx, nevt, 439 );
    storeWf( wfs, w_cx, nevt, 440 );
    storeWf( wfs, w_cx, nevt, 441 );
    storeWf( wfs, w_cx, nevt, 442 );
    storeWf( wfs, w_cx, nevt, 443 );
    storeWf( wfs, w_cx, nevt, 444 );
    storeWf( wfs, w_cx, nevt, 445 );
    storeWf( wfs, w_cx, nevt, 446 );
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
    storeWf( wfs, w_cx, nevt, 457 );
    storeWf( wfs, w_cx, nevt, 458 );
    storeWf( wfs, w_cx, nevt, 459 );
    storeWf( wfs, w_cx, nevt, 460 );
    storeWf( wfs, w_cx, nevt, 461 );
    storeWf( wfs, w_cx, nevt, 462 );
    storeWf( wfs, w_cx, nevt, 463 );
    storeWf( wfs, w_cx, nevt, 464 );
    storeWf( wfs, w_cx, nevt, 465 );
    storeWf( wfs, w_cx, nevt, 466 );
    storeWf( wfs, w_cx, nevt, 467 );
    storeWf( wfs, w_cx, nevt, 468 );
    storeWf( wfs, w_cx, nevt, 469 );
    storeWf( wfs, w_cx, nevt, 470 );
    storeWf( wfs, w_cx, nevt, 471 );
    storeWf( wfs, w_cx, nevt, 472 );
    storeWf( wfs, w_cx, nevt, 473 );
    storeWf( wfs, w_cx, nevt, 474 );
    storeWf( wfs, w_cx, nevt, 475 );
    storeWf( wfs, w_cx, nevt, 476 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
