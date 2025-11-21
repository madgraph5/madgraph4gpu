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
  diagramgroup5( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 8 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 15 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 30 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
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
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 72 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 74 );
    retrieveWf( wfs, w_cx, nevt, 78 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 80 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 83 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 89 );
    retrieveWf( wfs, w_cx, nevt, 91 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 93 );
    retrieveWf( wfs, w_cx, nevt, 94 );
    retrieveWf( wfs, w_cx, nevt, 95 );
    retrieveWf( wfs, w_cx, nevt, 96 );
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
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 128 );
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
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 155 );
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
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 229 );
    retrieveWf( wfs, w_cx, nevt, 236 );
    retrieveWf( wfs, w_cx, nevt, 237 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 240 );
    retrieveWf( wfs, w_cx, nevt, 241 );
    retrieveWf( wfs, w_cx, nevt, 242 );
    retrieveWf( wfs, w_cx, nevt, 244 );
#endif
#endif

    // *** DIAGRAM 801 OF 15495 ***
    // Wavefunction(s) for diagram number 801
    // (none)
    // Amplitude(s) for diagram number 801
    VVVV1_0( w_fp[37], w_fp[154], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    VVVV3_0( w_fp[37], w_fp[154], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVVV4_0( w_fp[37], w_fp[154], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 802 OF 15495 ***
    // Wavefunction(s) for diagram number 802
    // (none)
    // Amplitude(s) for diagram number 802
    VVV1_0( w_fp[154], w_fp[7], w_fp[40], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 803 OF 15495 ***
    // Wavefunction(s) for diagram number 803
    // (none)
    // Amplitude(s) for diagram number 803
    VVV1_0( w_fp[154], w_fp[4], w_fp[39], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 804 OF 15495 ***
    // Wavefunction(s) for diagram number 804
    // (none)
    // Amplitude(s) for diagram number 804
    FFV1_0( w_fp[218], w_fp[241], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];

    // *** DIAGRAM 805 OF 15495 ***
    // Wavefunction(s) for diagram number 805
    // (none)
    // Amplitude(s) for diagram number 805
    FFV1_0( w_fp[218], w_fp[2], w_fp[39], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 806 OF 15495 ***
    // Wavefunction(s) for diagram number 806
    // (none)
    // Amplitude(s) for diagram number 806
    FFV1_0( w_fp[171], w_fp[241], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];

    // *** DIAGRAM 807 OF 15495 ***
    // Wavefunction(s) for diagram number 807
    // (none)
    // Amplitude(s) for diagram number 807
    FFV1_0( w_fp[171], w_fp[2], w_fp[40], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 808 OF 15495 ***
    // Wavefunction(s) for diagram number 808
    // (none)
    // Amplitude(s) for diagram number 808
    VVVV1_0( w_fp[44], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    VVVV3_0( w_fp[44], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    VVVV4_0( w_fp[44], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 809 OF 15495 ***
    // Wavefunction(s) for diagram number 809
    // (none)
    // Amplitude(s) for diagram number 809
    VVV1_0( w_fp[154], w_fp[6], w_fp[47], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];

    // *** DIAGRAM 810 OF 15495 ***
    // Wavefunction(s) for diagram number 810
    // (none)
    // Amplitude(s) for diagram number 810
    VVV1_0( w_fp[154], w_fp[4], w_fp[46], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 811 OF 15495 ***
    // Wavefunction(s) for diagram number 811
    // (none)
    // Amplitude(s) for diagram number 811
    FFV1_0( w_fp[218], w_fp[242], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 812 OF 15495 ***
    // Wavefunction(s) for diagram number 812
    // (none)
    // Amplitude(s) for diagram number 812
    FFV1_0( w_fp[218], w_fp[2], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 813 OF 15495 ***
    // Wavefunction(s) for diagram number 813
    // (none)
    // Amplitude(s) for diagram number 813
    FFV1_0( w_fp[170], w_fp[242], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];

    // *** DIAGRAM 814 OF 15495 ***
    // Wavefunction(s) for diagram number 814
    // (none)
    // Amplitude(s) for diagram number 814
    FFV1_0( w_fp[170], w_fp[2], w_fp[47], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 815 OF 15495 ***
    // Wavefunction(s) for diagram number 815
    // (none)
    // Amplitude(s) for diagram number 815
    VVV1_0( w_fp[51], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVV1_0( w_fp[52], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    VVV1_0( w_fp[53], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 816 OF 15495 ***
    // Wavefunction(s) for diagram number 816
    // (none)
    // Amplitude(s) for diagram number 816
    FFV1_0( w_fp[171], w_fp[2], w_fp[51], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 817 OF 15495 ***
    // Wavefunction(s) for diagram number 817
    // (none)
    // Amplitude(s) for diagram number 817
    VVV1_0( w_fp[54], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVV1_0( w_fp[55], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    VVV1_0( w_fp[56], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 818 OF 15495 ***
    // Wavefunction(s) for diagram number 818
    // (none)
    // Amplitude(s) for diagram number 818
    FFV1_0( w_fp[170], w_fp[2], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 819 OF 15495 ***
    // Wavefunction(s) for diagram number 819
    // (none)
    // Amplitude(s) for diagram number 819
    VVV1_0( w_fp[63], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    VVV1_0( w_fp[64], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    VVV1_0( w_fp[65], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 820 OF 15495 ***
    // Wavefunction(s) for diagram number 820
    // (none)
    // Amplitude(s) for diagram number 820
    FFV1_0( w_fp[218], w_fp[2], w_fp[63], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 821 OF 15495 ***
    // Wavefunction(s) for diagram number 821
    FFV1_2( w_fp[168], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[200] );
    // Amplitude(s) for diagram number 821
    FFV1_0( w_fp[200], w_fp[229], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 822 OF 15495 ***
    // Wavefunction(s) for diagram number 822
    // (none)
    // Amplitude(s) for diagram number 822
    FFV1_0( w_fp[171], w_fp[229], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];

    // *** DIAGRAM 823 OF 15495 ***
    // Wavefunction(s) for diagram number 823
    // (none)
    // Amplitude(s) for diagram number 823
    FFV1_0( w_fp[168], w_fp[229], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 824 OF 15495 ***
    // Wavefunction(s) for diagram number 824
    FFV1_1( w_fp[2], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[118] );
    // Amplitude(s) for diagram number 824
    FFV1_0( w_fp[172], w_fp[118], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];

    // *** DIAGRAM 825 OF 15495 ***
    // Wavefunction(s) for diagram number 825
    // (none)
    // Amplitude(s) for diagram number 825
    FFV1_0( w_fp[172], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 826 OF 15495 ***
    // Wavefunction(s) for diagram number 826
    // (none)
    // Amplitude(s) for diagram number 826
    VVV1_0( w_fp[89], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 827 OF 15495 ***
    // Wavefunction(s) for diagram number 827
    // (none)
    // Amplitude(s) for diagram number 827
    FFV1_0( w_fp[171], w_fp[2], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 828 OF 15495 ***
    // Wavefunction(s) for diagram number 828
    // (none)
    // Amplitude(s) for diagram number 828
    VVV1_0( w_fp[44], w_fp[154], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 829 OF 15495 ***
    // Wavefunction(s) for diagram number 829
    // (none)
    // Amplitude(s) for diagram number 829
    FFV1_0( w_fp[168], w_fp[118], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 830 OF 15495 ***
    // Wavefunction(s) for diagram number 830
    // (none)
    // Amplitude(s) for diagram number 830
    FFV1_0( w_fp[200], w_fp[2], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 831 OF 15495 ***
    // Wavefunction(s) for diagram number 831
    // (none)
    // Amplitude(s) for diagram number 831
    VVV1_0( w_fp[8], w_fp[154], w_fp[88], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 832 OF 15495 ***
    // Wavefunction(s) for diagram number 832
    // (none)
    // Amplitude(s) for diagram number 832
    FFV1_0( w_fp[171], w_fp[118], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 833 OF 15495 ***
    // Wavefunction(s) for diagram number 833
    // (none)
    // Amplitude(s) for diagram number 833
    FFV1_0( w_fp[168], w_fp[2], w_fp[97], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[98], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[99], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 834 OF 15495 ***
    // Wavefunction(s) for diagram number 834
    FFV1_2( w_fp[168], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[99] );
    // Amplitude(s) for diagram number 834
    FFV1_0( w_fp[99], w_fp[229], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];

    // *** DIAGRAM 835 OF 15495 ***
    // Wavefunction(s) for diagram number 835
    // (none)
    // Amplitude(s) for diagram number 835
    FFV1_0( w_fp[170], w_fp[229], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];

    // *** DIAGRAM 836 OF 15495 ***
    // Wavefunction(s) for diagram number 836
    // (none)
    // Amplitude(s) for diagram number 836
    FFV1_0( w_fp[168], w_fp[229], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 837 OF 15495 ***
    // Wavefunction(s) for diagram number 837
    FFV1_1( w_fp[2], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[98] );
    // Amplitude(s) for diagram number 837
    FFV1_0( w_fp[172], w_fp[98], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 838 OF 15495 ***
    // Wavefunction(s) for diagram number 838
    // (none)
    // Amplitude(s) for diagram number 838
    FFV1_0( w_fp[172], w_fp[2], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 839 OF 15495 ***
    // Wavefunction(s) for diagram number 839
    // (none)
    // Amplitude(s) for diagram number 839
    VVV1_0( w_fp[105], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 840 OF 15495 ***
    // Wavefunction(s) for diagram number 840
    // (none)
    // Amplitude(s) for diagram number 840
    FFV1_0( w_fp[170], w_fp[2], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 841 OF 15495 ***
    // Wavefunction(s) for diagram number 841
    // (none)
    // Amplitude(s) for diagram number 841
    VVV1_0( w_fp[37], w_fp[154], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 842 OF 15495 ***
    // Wavefunction(s) for diagram number 842
    // (none)
    // Amplitude(s) for diagram number 842
    FFV1_0( w_fp[168], w_fp[98], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 843 OF 15495 ***
    // Wavefunction(s) for diagram number 843
    // (none)
    // Amplitude(s) for diagram number 843
    FFV1_0( w_fp[99], w_fp[2], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 844 OF 15495 ***
    // Wavefunction(s) for diagram number 844
    // (none)
    // Amplitude(s) for diagram number 844
    VVV1_0( w_fp[8], w_fp[154], w_fp[104], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];

    // *** DIAGRAM 845 OF 15495 ***
    // Wavefunction(s) for diagram number 845
    // (none)
    // Amplitude(s) for diagram number 845
    FFV1_0( w_fp[170], w_fp[98], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 846 OF 15495 ***
    // Wavefunction(s) for diagram number 846
    // (none)
    // Amplitude(s) for diagram number 846
    FFV1_0( w_fp[168], w_fp[2], w_fp[110], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[111], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[112], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 847 OF 15495 ***
    // Wavefunction(s) for diagram number 847
    // (none)
    // Amplitude(s) for diagram number 847
    FFV1_0( w_fp[218], w_fp[229], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 848 OF 15495 ***
    // Wavefunction(s) for diagram number 848
    FFV1_2( w_fp[168], w_fp[84], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[112] );
    // Amplitude(s) for diagram number 848
    FFV1_0( w_fp[112], w_fp[229], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];

    // *** DIAGRAM 849 OF 15495 ***
    // Wavefunction(s) for diagram number 849
    // (none)
    // Amplitude(s) for diagram number 849
    FFV1_0( w_fp[168], w_fp[229], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 850 OF 15495 ***
    // Wavefunction(s) for diagram number 850
    // (none)
    // Amplitude(s) for diagram number 850
    FFV1_0( w_fp[172], w_fp[128], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 851 OF 15495 ***
    // Wavefunction(s) for diagram number 851
    // (none)
    // Amplitude(s) for diagram number 851
    FFV1_0( w_fp[172], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 852 OF 15495 ***
    // Wavefunction(s) for diagram number 852
    // (none)
    // Amplitude(s) for diagram number 852
    VVV1_0( w_fp[10], w_fp[154], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 853 OF 15495 ***
    // Wavefunction(s) for diagram number 853
    // (none)
    // Amplitude(s) for diagram number 853
    FFV1_0( w_fp[168], w_fp[128], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 854 OF 15495 ***
    // Wavefunction(s) for diagram number 854
    // (none)
    // Amplitude(s) for diagram number 854
    FFV1_0( w_fp[112], w_fp[2], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 855 OF 15495 ***
    // Wavefunction(s) for diagram number 855
    // (none)
    // Amplitude(s) for diagram number 855
    VVV1_0( w_fp[85], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];

    // *** DIAGRAM 856 OF 15495 ***
    // Wavefunction(s) for diagram number 856
    // (none)
    // Amplitude(s) for diagram number 856
    FFV1_0( w_fp[218], w_fp[2], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 857 OF 15495 ***
    // Wavefunction(s) for diagram number 857
    // (none)
    // Amplitude(s) for diagram number 857
    VVV1_0( w_fp[8], w_fp[154], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 858 OF 15495 ***
    // Wavefunction(s) for diagram number 858
    // (none)
    // Amplitude(s) for diagram number 858
    FFV1_0( w_fp[218], w_fp[128], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];

    // *** DIAGRAM 859 OF 15495 ***
    // Wavefunction(s) for diagram number 859
    // (none)
    // Amplitude(s) for diagram number 859
    FFV1_0( w_fp[168], w_fp[2], w_fp[93], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[92], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[91], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 860 OF 15495 ***
    // Wavefunction(s) for diagram number 860
    // (none)
    // Amplitude(s) for diagram number 860
    FFV1_0( w_fp[168], w_fp[229], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[229], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[229], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 861 OF 15495 ***
    // Wavefunction(s) for diagram number 861
    // (none)
    // Amplitude(s) for diagram number 861
    FFV1_0( w_fp[172], w_fp[2], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[172], w_fp[2], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[172], w_fp[2], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 862 OF 15495 ***
    // Wavefunction(s) for diagram number 862
    // (none)
    // Amplitude(s) for diagram number 862
    FFV1_0( w_fp[168], w_fp[2], w_fp[148], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[149], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[150], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 863 OF 15495 ***
    // Wavefunction(s) for diagram number 863
    // (none)
    // Amplitude(s) for diagram number 863
    FFV1_0( w_fp[202], w_fp[238], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 864 OF 15495 ***
    // Wavefunction(s) for diagram number 864
    // (none)
    // Amplitude(s) for diagram number 864
    FFV1_0( w_fp[202], w_fp[237], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 865 OF 15495 ***
    // Wavefunction(s) for diagram number 865
    // (none)
    // Amplitude(s) for diagram number 865
    FFV1_0( w_fp[175], w_fp[155], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 866 OF 15495 ***
    // Wavefunction(s) for diagram number 866
    // (none)
    // Amplitude(s) for diagram number 866
    FFV1_0( w_fp[175], w_fp[237], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 867 OF 15495 ***
    // Wavefunction(s) for diagram number 867
    // (none)
    // Amplitude(s) for diagram number 867
    FFV1_0( w_fp[176], w_fp[155], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 868 OF 15495 ***
    // Wavefunction(s) for diagram number 868
    // (none)
    // Amplitude(s) for diagram number 868
    FFV1_0( w_fp[176], w_fp[238], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 869 OF 15495 ***
    // Wavefunction(s) for diagram number 869
    FFV1P0_3( w_fp[174], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[150] );
    // Amplitude(s) for diagram number 869
    VVVV1_0( w_fp[10], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV3_0( w_fp[10], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0( w_fp[10], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 870 OF 15495 ***
    // Wavefunction(s) for diagram number 870
    // (none)
    // Amplitude(s) for diagram number 870
    VVV1_0( w_fp[150], w_fp[7], w_fp[15], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 871 OF 15495 ***
    // Wavefunction(s) for diagram number 871
    // (none)
    // Amplitude(s) for diagram number 871
    VVV1_0( w_fp[150], w_fp[5], w_fp[13], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 872 OF 15495 ***
    // Wavefunction(s) for diagram number 872
    // (none)
    // Amplitude(s) for diagram number 872
    FFV1_0( w_fp[175], w_fp[67], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];

    // *** DIAGRAM 873 OF 15495 ***
    // Wavefunction(s) for diagram number 873
    // (none)
    // Amplitude(s) for diagram number 873
    FFV1_0( w_fp[175], w_fp[2], w_fp[13], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 874 OF 15495 ***
    // Wavefunction(s) for diagram number 874
    // (none)
    // Amplitude(s) for diagram number 874
    FFV1_0( w_fp[176], w_fp[67], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];

    // *** DIAGRAM 875 OF 15495 ***
    // Wavefunction(s) for diagram number 875
    // (none)
    // Amplitude(s) for diagram number 875
    FFV1_0( w_fp[176], w_fp[2], w_fp[15], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 876 OF 15495 ***
    // Wavefunction(s) for diagram number 876
    // (none)
    // Amplitude(s) for diagram number 876
    VVVV1_0( w_fp[26], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    VVVV3_0( w_fp[26], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVVV4_0( w_fp[26], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 877 OF 15495 ***
    // Wavefunction(s) for diagram number 877
    // (none)
    // Amplitude(s) for diagram number 877
    VVV1_0( w_fp[150], w_fp[7], w_fp[30], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 878 OF 15495 ***
    // Wavefunction(s) for diagram number 878
    // (none)
    // Amplitude(s) for diagram number 878
    VVV1_0( w_fp[150], w_fp[4], w_fp[29], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 879 OF 15495 ***
    // Wavefunction(s) for diagram number 879
    // (none)
    // Amplitude(s) for diagram number 879
    FFV1_0( w_fp[202], w_fp[240], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];

    // *** DIAGRAM 880 OF 15495 ***
    // Wavefunction(s) for diagram number 880
    // (none)
    // Amplitude(s) for diagram number 880
    FFV1_0( w_fp[202], w_fp[2], w_fp[29], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 881 OF 15495 ***
    // Wavefunction(s) for diagram number 881
    // (none)
    // Amplitude(s) for diagram number 881
    FFV1_0( w_fp[176], w_fp[240], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];

    // *** DIAGRAM 882 OF 15495 ***
    // Wavefunction(s) for diagram number 882
    // (none)
    // Amplitude(s) for diagram number 882
    FFV1_0( w_fp[176], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 883 OF 15495 ***
    // Wavefunction(s) for diagram number 883
    // (none)
    // Amplitude(s) for diagram number 883
    VVVV1_0( w_fp[44], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[44], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    VVVV4_0( w_fp[44], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 884 OF 15495 ***
    // Wavefunction(s) for diagram number 884
    // (none)
    // Amplitude(s) for diagram number 884
    VVV1_0( w_fp[150], w_fp[5], w_fp[47], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];

    // *** DIAGRAM 885 OF 15495 ***
    // Wavefunction(s) for diagram number 885
    // (none)
    // Amplitude(s) for diagram number 885
    VVV1_0( w_fp[150], w_fp[4], w_fp[45], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 886 OF 15495 ***
    // Wavefunction(s) for diagram number 886
    // (none)
    // Amplitude(s) for diagram number 886
    FFV1_0( w_fp[202], w_fp[242], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 887 OF 15495 ***
    // Wavefunction(s) for diagram number 887
    // (none)
    // Amplitude(s) for diagram number 887
    FFV1_0( w_fp[202], w_fp[2], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 888 OF 15495 ***
    // Wavefunction(s) for diagram number 888
    // (none)
    // Amplitude(s) for diagram number 888
    FFV1_0( w_fp[175], w_fp[242], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];

    // *** DIAGRAM 889 OF 15495 ***
    // Wavefunction(s) for diagram number 889
    // (none)
    // Amplitude(s) for diagram number 889
    FFV1_0( w_fp[175], w_fp[2], w_fp[47], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 890 OF 15495 ***
    // Wavefunction(s) for diagram number 890
    // (none)
    // Amplitude(s) for diagram number 890
    VVV1_0( w_fp[48], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVV1_0( w_fp[49], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    VVV1_0( w_fp[50], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 891 OF 15495 ***
    // Wavefunction(s) for diagram number 891
    // (none)
    // Amplitude(s) for diagram number 891
    FFV1_0( w_fp[176], w_fp[2], w_fp[48], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[49], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[50], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 892 OF 15495 ***
    // Wavefunction(s) for diagram number 892
    // (none)
    // Amplitude(s) for diagram number 892
    VVV1_0( w_fp[54], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    VVV1_0( w_fp[55], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    VVV1_0( w_fp[56], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 893 OF 15495 ***
    // Wavefunction(s) for diagram number 893
    // (none)
    // Amplitude(s) for diagram number 893
    FFV1_0( w_fp[175], w_fp[2], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 894 OF 15495 ***
    // Wavefunction(s) for diagram number 894
    // (none)
    // Amplitude(s) for diagram number 894
    VVV1_0( w_fp[60], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVV1_0( w_fp[61], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    VVV1_0( w_fp[62], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 895 OF 15495 ***
    // Wavefunction(s) for diagram number 895
    // (none)
    // Amplitude(s) for diagram number 895
    FFV1_0( w_fp[202], w_fp[2], w_fp[60], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[202], w_fp[2], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[202], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 896 OF 15495 ***
    // Wavefunction(s) for diagram number 896
    FFV1_2( w_fp[174], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[149] );
    // Amplitude(s) for diagram number 896
    FFV1_0( w_fp[149], w_fp[229], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];

    // *** DIAGRAM 897 OF 15495 ***
    // Wavefunction(s) for diagram number 897
    // (none)
    // Amplitude(s) for diagram number 897
    FFV1_0( w_fp[176], w_fp[229], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];

    // *** DIAGRAM 898 OF 15495 ***
    // Wavefunction(s) for diagram number 898
    // (none)
    // Amplitude(s) for diagram number 898
    FFV1_0( w_fp[174], w_fp[229], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 899 OF 15495 ***
    // Wavefunction(s) for diagram number 899
    FFV1_1( w_fp[2], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[148] );
    // Amplitude(s) for diagram number 899
    FFV1_0( w_fp[177], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];

    // *** DIAGRAM 900 OF 15495 ***
    // Wavefunction(s) for diagram number 900
    // (none)
    // Amplitude(s) for diagram number 900
    FFV1_0( w_fp[177], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 901 OF 15495 ***
    // Wavefunction(s) for diagram number 901
    // (none)
    // Amplitude(s) for diagram number 901
    VVV1_0( w_fp[70], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 902 OF 15495 ***
    // Wavefunction(s) for diagram number 902
    // (none)
    // Amplitude(s) for diagram number 902
    FFV1_0( w_fp[176], w_fp[2], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 903 OF 15495 ***
    // Wavefunction(s) for diagram number 903
    // (none)
    // Amplitude(s) for diagram number 903
    VVV1_0( w_fp[44], w_fp[150], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 904 OF 15495 ***
    // Wavefunction(s) for diagram number 904
    // (none)
    // Amplitude(s) for diagram number 904
    FFV1_0( w_fp[174], w_fp[148], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 905 OF 15495 ***
    // Wavefunction(s) for diagram number 905
    // (none)
    // Amplitude(s) for diagram number 905
    FFV1_0( w_fp[149], w_fp[2], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 906 OF 15495 ***
    // Wavefunction(s) for diagram number 906
    // (none)
    // Amplitude(s) for diagram number 906
    VVV1_0( w_fp[8], w_fp[150], w_fp[69], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 907 OF 15495 ***
    // Wavefunction(s) for diagram number 907
    // (none)
    // Amplitude(s) for diagram number 907
    FFV1_0( w_fp[176], w_fp[148], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

    // *** DIAGRAM 908 OF 15495 ***
    // Wavefunction(s) for diagram number 908
    // (none)
    // Amplitude(s) for diagram number 908
    FFV1_0( w_fp[174], w_fp[2], w_fp[81], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[82], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[174], w_fp[2], w_fp[83], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 909 OF 15495 ***
    // Wavefunction(s) for diagram number 909
    FFV1_2( w_fp[174], w_fp[102], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[83] );
    // Amplitude(s) for diagram number 909
    FFV1_0( w_fp[83], w_fp[229], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];

    // *** DIAGRAM 910 OF 15495 ***
    // Wavefunction(s) for diagram number 910
    // (none)
    // Amplitude(s) for diagram number 910
    FFV1_0( w_fp[175], w_fp[229], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];

    // *** DIAGRAM 911 OF 15495 ***
    // Wavefunction(s) for diagram number 911
    // (none)
    // Amplitude(s) for diagram number 911
    FFV1_0( w_fp[174], w_fp[229], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 912 OF 15495 ***
    // Wavefunction(s) for diagram number 912
    // (none)
    // Amplitude(s) for diagram number 912
    FFV1_0( w_fp[177], w_fp[98], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];

    // *** DIAGRAM 913 OF 15495 ***
    // Wavefunction(s) for diagram number 913
    // (none)
    // Amplitude(s) for diagram number 913
    FFV1_0( w_fp[177], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 914 OF 15495 ***
    // Wavefunction(s) for diagram number 914
    // (none)
    // Amplitude(s) for diagram number 914
    VVV1_0( w_fp[105], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 915 OF 15495 ***
    // Wavefunction(s) for diagram number 915
    // (none)
    // Amplitude(s) for diagram number 915
    FFV1_0( w_fp[175], w_fp[2], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 916 OF 15495 ***
    // Wavefunction(s) for diagram number 916
    // (none)
    // Amplitude(s) for diagram number 916
    VVV1_0( w_fp[26], w_fp[150], w_fp[102], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];

    // *** DIAGRAM 917 OF 15495 ***
    // Wavefunction(s) for diagram number 917
    // (none)
    // Amplitude(s) for diagram number 917
    FFV1_0( w_fp[174], w_fp[98], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 918 OF 15495 ***
    // Wavefunction(s) for diagram number 918
    // (none)
    // Amplitude(s) for diagram number 918
    FFV1_0( w_fp[83], w_fp[2], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 919 OF 15495 ***
    // Wavefunction(s) for diagram number 919
    // (none)
    // Amplitude(s) for diagram number 919
    VVV1_0( w_fp[8], w_fp[150], w_fp[103], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];

    // *** DIAGRAM 920 OF 15495 ***
    // Wavefunction(s) for diagram number 920
    // (none)
    // Amplitude(s) for diagram number 920
    FFV1_0( w_fp[175], w_fp[98], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 921 OF 15495 ***
    // Wavefunction(s) for diagram number 921
    // (none)
    // Amplitude(s) for diagram number 921
    FFV1_0( w_fp[174], w_fp[2], w_fp[107], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[108], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[109], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 922 OF 15495 ***
    // Wavefunction(s) for diagram number 922
    // (none)
    // Amplitude(s) for diagram number 922
    FFV1_0( w_fp[202], w_fp[229], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];

    // *** DIAGRAM 923 OF 15495 ***
    // Wavefunction(s) for diagram number 923
    FFV1_2( w_fp[174], w_fp[100], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[109] );
    // Amplitude(s) for diagram number 923
    FFV1_0( w_fp[109], w_fp[229], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];

    // *** DIAGRAM 924 OF 15495 ***
    // Wavefunction(s) for diagram number 924
    // (none)
    // Amplitude(s) for diagram number 924
    FFV1_0( w_fp[174], w_fp[229], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 925 OF 15495 ***
    // Wavefunction(s) for diagram number 925
    // (none)
    // Amplitude(s) for diagram number 925
    FFV1_0( w_fp[177], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 926 OF 15495 ***
    // Wavefunction(s) for diagram number 926
    // (none)
    // Amplitude(s) for diagram number 926
    FFV1_0( w_fp[177], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 927 OF 15495 ***
    // Wavefunction(s) for diagram number 927
    // (none)
    // Amplitude(s) for diagram number 927
    VVV1_0( w_fp[10], w_fp[150], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 928 OF 15495 ***
    // Wavefunction(s) for diagram number 928
    // (none)
    // Amplitude(s) for diagram number 928
    FFV1_0( w_fp[174], w_fp[122], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 929 OF 15495 ***
    // Wavefunction(s) for diagram number 929
    // (none)
    // Amplitude(s) for diagram number 929
    FFV1_0( w_fp[109], w_fp[2], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 930 OF 15495 ***
    // Wavefunction(s) for diagram number 930
    // (none)
    // Amplitude(s) for diagram number 930
    VVV1_0( w_fp[101], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 931 OF 15495 ***
    // Wavefunction(s) for diagram number 931
    // (none)
    // Amplitude(s) for diagram number 931
    FFV1_0( w_fp[202], w_fp[2], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 932 OF 15495 ***
    // Wavefunction(s) for diagram number 932
    // (none)
    // Amplitude(s) for diagram number 932
    VVV1_0( w_fp[8], w_fp[150], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 933 OF 15495 ***
    // Wavefunction(s) for diagram number 933
    // (none)
    // Amplitude(s) for diagram number 933
    FFV1_0( w_fp[202], w_fp[122], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 934 OF 15495 ***
    // Wavefunction(s) for diagram number 934
    // (none)
    // Amplitude(s) for diagram number 934
    FFV1_0( w_fp[174], w_fp[2], w_fp[74], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[73], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[72], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 935 OF 15495 ***
    // Wavefunction(s) for diagram number 935
    // (none)
    // Amplitude(s) for diagram number 935
    FFV1_0( w_fp[174], w_fp[229], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[229], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[229], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 936 OF 15495 ***
    // Wavefunction(s) for diagram number 936
    // (none)
    // Amplitude(s) for diagram number 936
    FFV1_0( w_fp[177], w_fp[2], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[177], w_fp[2], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[177], w_fp[2], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 937 OF 15495 ***
    // Wavefunction(s) for diagram number 937
    // (none)
    // Amplitude(s) for diagram number 937
    FFV1_0( w_fp[174], w_fp[2], w_fp[142], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[143], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[144], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 938 OF 15495 ***
    // Wavefunction(s) for diagram number 938
    // (none)
    // Amplitude(s) for diagram number 938
    FFV1_0( w_fp[204], w_fp[238], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 939 OF 15495 ***
    // Wavefunction(s) for diagram number 939
    // (none)
    // Amplitude(s) for diagram number 939
    FFV1_0( w_fp[204], w_fp[236], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 940 OF 15495 ***
    // Wavefunction(s) for diagram number 940
    // (none)
    // Amplitude(s) for diagram number 940
    FFV1_0( w_fp[180], w_fp[155], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 941 OF 15495 ***
    // Wavefunction(s) for diagram number 941
    // (none)
    // Amplitude(s) for diagram number 941
    FFV1_0( w_fp[180], w_fp[236], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 942 OF 15495 ***
    // Wavefunction(s) for diagram number 942
    // (none)
    // Amplitude(s) for diagram number 942
    FFV1_0( w_fp[181], w_fp[155], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 943 OF 15495 ***
    // Wavefunction(s) for diagram number 943
    // (none)
    // Amplitude(s) for diagram number 943
    FFV1_0( w_fp[181], w_fp[238], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 944 OF 15495 ***
    // Wavefunction(s) for diagram number 944
    FFV1P0_3( w_fp[179], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[144] );
    // Amplitude(s) for diagram number 944
    VVVV1_0( w_fp[10], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV3_0( w_fp[10], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0( w_fp[10], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 945 OF 15495 ***
    // Wavefunction(s) for diagram number 945
    // (none)
    // Amplitude(s) for diagram number 945
    VVV1_0( w_fp[144], w_fp[6], w_fp[15], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 946 OF 15495 ***
    // Wavefunction(s) for diagram number 946
    // (none)
    // Amplitude(s) for diagram number 946
    VVV1_0( w_fp[144], w_fp[5], w_fp[12], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 947 OF 15495 ***
    // Wavefunction(s) for diagram number 947
    // (none)
    // Amplitude(s) for diagram number 947
    FFV1_0( w_fp[180], w_fp[67], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];

    // *** DIAGRAM 948 OF 15495 ***
    // Wavefunction(s) for diagram number 948
    // (none)
    // Amplitude(s) for diagram number 948
    FFV1_0( w_fp[180], w_fp[2], w_fp[12], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 949 OF 15495 ***
    // Wavefunction(s) for diagram number 949
    // (none)
    // Amplitude(s) for diagram number 949
    FFV1_0( w_fp[181], w_fp[67], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];

    // *** DIAGRAM 950 OF 15495 ***
    // Wavefunction(s) for diagram number 950
    // (none)
    // Amplitude(s) for diagram number 950
    FFV1_0( w_fp[181], w_fp[2], w_fp[15], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 951 OF 15495 ***
    // Wavefunction(s) for diagram number 951
    // (none)
    // Amplitude(s) for diagram number 951
    VVVV1_0( w_fp[26], w_fp[144], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    VVVV3_0( w_fp[26], w_fp[144], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    VVVV4_0( w_fp[26], w_fp[144], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 952 OF 15495 ***
    // Wavefunction(s) for diagram number 952
    // (none)
    // Amplitude(s) for diagram number 952
    VVV1_0( w_fp[144], w_fp[6], w_fp[30], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 953 OF 15495 ***
    // Wavefunction(s) for diagram number 953
    // (none)
    // Amplitude(s) for diagram number 953
    VVV1_0( w_fp[144], w_fp[4], w_fp[28], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 954 OF 15495 ***
    // Wavefunction(s) for diagram number 954
    // (none)
    // Amplitude(s) for diagram number 954
    FFV1_0( w_fp[204], w_fp[240], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];

    // *** DIAGRAM 955 OF 15495 ***
    // Wavefunction(s) for diagram number 955
    // (none)
    // Amplitude(s) for diagram number 955
    FFV1_0( w_fp[204], w_fp[2], w_fp[28], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 956 OF 15495 ***
    // Wavefunction(s) for diagram number 956
    // (none)
    // Amplitude(s) for diagram number 956
    FFV1_0( w_fp[181], w_fp[240], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];

    // *** DIAGRAM 957 OF 15495 ***
    // Wavefunction(s) for diagram number 957
    // (none)
    // Amplitude(s) for diagram number 957
    FFV1_0( w_fp[181], w_fp[2], w_fp[30], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 958 OF 15495 ***
    // Wavefunction(s) for diagram number 958
    // (none)
    // Amplitude(s) for diagram number 958
    VVVV1_0( w_fp[37], w_fp[144], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    VVVV3_0( w_fp[37], w_fp[144], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    VVVV4_0( w_fp[37], w_fp[144], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 959 OF 15495 ***
    // Wavefunction(s) for diagram number 959
    // (none)
    // Amplitude(s) for diagram number 959
    VVV1_0( w_fp[144], w_fp[5], w_fp[40], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];

    // *** DIAGRAM 960 OF 15495 ***
    // Wavefunction(s) for diagram number 960
    // (none)
    // Amplitude(s) for diagram number 960
    VVV1_0( w_fp[144], w_fp[4], w_fp[38], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 961 OF 15495 ***
    // Wavefunction(s) for diagram number 961
    // (none)
    // Amplitude(s) for diagram number 961
    FFV1_0( w_fp[204], w_fp[241], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];

    // *** DIAGRAM 962 OF 15495 ***
    // Wavefunction(s) for diagram number 962
    // (none)
    // Amplitude(s) for diagram number 962
    FFV1_0( w_fp[204], w_fp[2], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 963 OF 15495 ***
    // Wavefunction(s) for diagram number 963
    // (none)
    // Amplitude(s) for diagram number 963
    FFV1_0( w_fp[180], w_fp[241], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];

    // *** DIAGRAM 964 OF 15495 ***
    // Wavefunction(s) for diagram number 964
    // (none)
    // Amplitude(s) for diagram number 964
    FFV1_0( w_fp[180], w_fp[2], w_fp[40], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 965 OF 15495 ***
    // Wavefunction(s) for diagram number 965
    // (none)
    // Amplitude(s) for diagram number 965
    VVV1_0( w_fp[48], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVV1_0( w_fp[49], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    VVV1_0( w_fp[50], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];

    // *** DIAGRAM 966 OF 15495 ***
    // Wavefunction(s) for diagram number 966
    // (none)
    // Amplitude(s) for diagram number 966
    FFV1_0( w_fp[181], w_fp[2], w_fp[48], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[49], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[50], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 967 OF 15495 ***
    // Wavefunction(s) for diagram number 967
    // (none)
    // Amplitude(s) for diagram number 967
    VVV1_0( w_fp[51], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    VVV1_0( w_fp[52], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    VVV1_0( w_fp[53], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];

    // *** DIAGRAM 968 OF 15495 ***
    // Wavefunction(s) for diagram number 968
    // (none)
    // Amplitude(s) for diagram number 968
    FFV1_0( w_fp[180], w_fp[2], w_fp[51], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[180], w_fp[2], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[180], w_fp[2], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 969 OF 15495 ***
    // Wavefunction(s) for diagram number 969
    // (none)
    // Amplitude(s) for diagram number 969
    VVV1_0( w_fp[57], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    VVV1_0( w_fp[58], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    VVV1_0( w_fp[59], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];

    // *** DIAGRAM 970 OF 15495 ***
    // Wavefunction(s) for diagram number 970
    // (none)
    // Amplitude(s) for diagram number 970
    FFV1_0( w_fp[204], w_fp[2], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[204], w_fp[2], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[204], w_fp[2], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 971 OF 15495 ***
    // Wavefunction(s) for diagram number 971
    FFV1_2( w_fp[179], w_fp[66], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[143] );
    // Amplitude(s) for diagram number 971
    FFV1_0( w_fp[143], w_fp[229], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];

    // *** DIAGRAM 972 OF 15495 ***
    // Wavefunction(s) for diagram number 972
    // (none)
    // Amplitude(s) for diagram number 972
    FFV1_0( w_fp[181], w_fp[229], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];

    // *** DIAGRAM 973 OF 15495 ***
    // Wavefunction(s) for diagram number 973
    // (none)
    // Amplitude(s) for diagram number 973
    FFV1_0( w_fp[179], w_fp[229], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 974 OF 15495 ***
    // Wavefunction(s) for diagram number 974
    // (none)
    // Amplitude(s) for diagram number 974
    FFV1_0( w_fp[182], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];

    // *** DIAGRAM 975 OF 15495 ***
    // Wavefunction(s) for diagram number 975
    // (none)
    // Amplitude(s) for diagram number 975
    FFV1_0( w_fp[182], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 976 OF 15495 ***
    // Wavefunction(s) for diagram number 976
    // (none)
    // Amplitude(s) for diagram number 976
    VVV1_0( w_fp[70], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 977 OF 15495 ***
    // Wavefunction(s) for diagram number 977
    // (none)
    // Amplitude(s) for diagram number 977
    FFV1_0( w_fp[181], w_fp[2], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 978 OF 15495 ***
    // Wavefunction(s) for diagram number 978
    // (none)
    // Amplitude(s) for diagram number 978
    VVV1_0( w_fp[37], w_fp[144], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];

    // *** DIAGRAM 979 OF 15495 ***
    // Wavefunction(s) for diagram number 979
    // (none)
    // Amplitude(s) for diagram number 979
    FFV1_0( w_fp[179], w_fp[148], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 980 OF 15495 ***
    // Wavefunction(s) for diagram number 980
    // (none)
    // Amplitude(s) for diagram number 980
    FFV1_0( w_fp[143], w_fp[2], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 981 OF 15495 ***
    // Wavefunction(s) for diagram number 981
    // (none)
    // Amplitude(s) for diagram number 981
    VVV1_0( w_fp[8], w_fp[144], w_fp[68], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 982 OF 15495 ***
    // Wavefunction(s) for diagram number 982
    // (none)
    // Amplitude(s) for diagram number 982
    FFV1_0( w_fp[181], w_fp[148], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];

    // *** DIAGRAM 983 OF 15495 ***
    // Wavefunction(s) for diagram number 983
    // (none)
    // Amplitude(s) for diagram number 983
    FFV1_0( w_fp[179], w_fp[2], w_fp[78], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[79], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[80], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 984 OF 15495 ***
    // Wavefunction(s) for diagram number 984
    FFV1_2( w_fp[179], w_fp[86], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[80] );
    // Amplitude(s) for diagram number 984
    FFV1_0( w_fp[80], w_fp[229], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];

    // *** DIAGRAM 985 OF 15495 ***
    // Wavefunction(s) for diagram number 985
    // (none)
    // Amplitude(s) for diagram number 985
    FFV1_0( w_fp[180], w_fp[229], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];

    // *** DIAGRAM 986 OF 15495 ***
    // Wavefunction(s) for diagram number 986
    // (none)
    // Amplitude(s) for diagram number 986
    FFV1_0( w_fp[179], w_fp[229], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 987 OF 15495 ***
    // Wavefunction(s) for diagram number 987
    // (none)
    // Amplitude(s) for diagram number 987
    FFV1_0( w_fp[182], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];

    // *** DIAGRAM 988 OF 15495 ***
    // Wavefunction(s) for diagram number 988
    // (none)
    // Amplitude(s) for diagram number 988
    FFV1_0( w_fp[182], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 989 OF 15495 ***
    // Wavefunction(s) for diagram number 989
    // (none)
    // Amplitude(s) for diagram number 989
    VVV1_0( w_fp[89], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 990 OF 15495 ***
    // Wavefunction(s) for diagram number 990
    // (none)
    // Amplitude(s) for diagram number 990
    FFV1_0( w_fp[180], w_fp[2], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 991 OF 15495 ***
    // Wavefunction(s) for diagram number 991
    // (none)
    // Amplitude(s) for diagram number 991
    VVV1_0( w_fp[26], w_fp[144], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];

    // *** DIAGRAM 992 OF 15495 ***
    // Wavefunction(s) for diagram number 992
    // (none)
    // Amplitude(s) for diagram number 992
    FFV1_0( w_fp[179], w_fp[118], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 993 OF 15495 ***
    // Wavefunction(s) for diagram number 993
    // (none)
    // Amplitude(s) for diagram number 993
    FFV1_0( w_fp[80], w_fp[2], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 994 OF 15495 ***
    // Wavefunction(s) for diagram number 994
    // (none)
    // Amplitude(s) for diagram number 994
    VVV1_0( w_fp[8], w_fp[144], w_fp[87], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 995 OF 15495 ***
    // Wavefunction(s) for diagram number 995
    // (none)
    // Amplitude(s) for diagram number 995
    FFV1_0( w_fp[180], w_fp[118], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];

    // *** DIAGRAM 996 OF 15495 ***
    // Wavefunction(s) for diagram number 996
    // (none)
    // Amplitude(s) for diagram number 996
    FFV1_0( w_fp[179], w_fp[2], w_fp[94], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[95], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[96], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 997 OF 15495 ***
    // Wavefunction(s) for diagram number 997
    // (none)
    // Amplitude(s) for diagram number 997
    FFV1_0( w_fp[204], w_fp[229], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];

    // *** DIAGRAM 998 OF 15495 ***
    // Wavefunction(s) for diagram number 998
    FFV1_2( w_fp[179], w_fp[113], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[96] );
    // Amplitude(s) for diagram number 998
    FFV1_0( w_fp[96], w_fp[229], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];

    // *** DIAGRAM 999 OF 15495 ***
    // Wavefunction(s) for diagram number 999
    // (none)
    // Amplitude(s) for diagram number 999
    FFV1_0( w_fp[179], w_fp[229], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1000 OF 15495 ***
    // Wavefunction(s) for diagram number 1000
    // (none)
    // Amplitude(s) for diagram number 1000
    FFV1_0( w_fp[182], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 80 );
    storeWf( wfs, w_cx, nevt, 83 );
    storeWf( wfs, w_cx, nevt, 96 );
    storeWf( wfs, w_cx, nevt, 98 );
    storeWf( wfs, w_cx, nevt, 99 );
    storeWf( wfs, w_cx, nevt, 109 );
    storeWf( wfs, w_cx, nevt, 112 );
    storeWf( wfs, w_cx, nevt, 118 );
    storeWf( wfs, w_cx, nevt, 143 );
    storeWf( wfs, w_cx, nevt, 144 );
    storeWf( wfs, w_cx, nevt, 148 );
    storeWf( wfs, w_cx, nevt, 149 );
    storeWf( wfs, w_cx, nevt, 150 );
    storeWf( wfs, w_cx, nevt, 200 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
