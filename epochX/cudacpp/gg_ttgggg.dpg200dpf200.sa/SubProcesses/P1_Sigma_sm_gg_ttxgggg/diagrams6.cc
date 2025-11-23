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
  diagramgroup6( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 47 );
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
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 89 );
    retrieveWf( wfs, w_cx, nevt, 96 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 105 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 138 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 155 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 165 );
    retrieveWf( wfs, w_cx, nevt, 166 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 192 );
    retrieveWf( wfs, w_cx, nevt, 195 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 221 );
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

    // *** DIAGRAM 1001 OF 15495 ***
    // Wavefunction(s) for diagram number 1001
    // (none)
    // Amplitude(s) for diagram number 1001
    FFV1_0( w_fp[182], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1002 OF 15495 ***
    // Wavefunction(s) for diagram number 1002
    // (none)
    // Amplitude(s) for diagram number 1002
    VVV1_0( w_fp[10], w_fp[144], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1003 OF 15495 ***
    // Wavefunction(s) for diagram number 1003
    // (none)
    // Amplitude(s) for diagram number 1003
    FFV1_0( w_fp[179], w_fp[244], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1004 OF 15495 ***
    // Wavefunction(s) for diagram number 1004
    // (none)
    // Amplitude(s) for diagram number 1004
    FFV1_0( w_fp[96], w_fp[2], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1005 OF 15495 ***
    // Wavefunction(s) for diagram number 1005
    // (none)
    // Amplitude(s) for diagram number 1005
    VVV1_0( w_fp[114], w_fp[144], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1006 OF 15495 ***
    // Wavefunction(s) for diagram number 1006
    // (none)
    // Amplitude(s) for diagram number 1006
    FFV1_0( w_fp[204], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1007 OF 15495 ***
    // Wavefunction(s) for diagram number 1007
    // (none)
    // Amplitude(s) for diagram number 1007
    VVV1_0( w_fp[8], w_fp[144], w_fp[115], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 1008 OF 15495 ***
    // Wavefunction(s) for diagram number 1008
    // (none)
    // Amplitude(s) for diagram number 1008
    FFV1_0( w_fp[204], w_fp[244], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];

    // *** DIAGRAM 1009 OF 15495 ***
    // Wavefunction(s) for diagram number 1009
    // (none)
    // Amplitude(s) for diagram number 1009
    FFV1_0( w_fp[179], w_fp[2], w_fp[77], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[76], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[179], w_fp[2], w_fp[75], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1010 OF 15495 ***
    // Wavefunction(s) for diagram number 1010
    // (none)
    // Amplitude(s) for diagram number 1010
    FFV1_0( w_fp[179], w_fp[229], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[179], w_fp[229], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[229], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1011 OF 15495 ***
    // Wavefunction(s) for diagram number 1011
    // (none)
    // Amplitude(s) for diagram number 1011
    FFV1_0( w_fp[182], w_fp[2], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[182], w_fp[2], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[182], w_fp[2], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1012 OF 15495 ***
    // Wavefunction(s) for diagram number 1012
    // (none)
    // Amplitude(s) for diagram number 1012
    FFV1_0( w_fp[179], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[179], w_fp[2], w_fp[138], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1013 OF 15495 ***
    // Wavefunction(s) for diagram number 1013
    // (none)
    // Amplitude(s) for diagram number 1013
    FFV1_0( w_fp[221], w_fp[236], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];

    // *** DIAGRAM 1014 OF 15495 ***
    // Wavefunction(s) for diagram number 1014
    // (none)
    // Amplitude(s) for diagram number 1014
    FFV1_0( w_fp[221], w_fp[237], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 1015 OF 15495 ***
    // Wavefunction(s) for diagram number 1015
    FFV1P0_3( w_fp[3], w_fp[229], COUPs[1], 1.0, depCoup, 0., 0., w_fp[138] );
    // Amplitude(s) for diagram number 1015
    VVV1_0( w_fp[68], w_fp[7], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 1016 OF 15495 ***
    // Wavefunction(s) for diagram number 1016
    // (none)
    // Amplitude(s) for diagram number 1016
    FFV1_0( w_fp[3], w_fp[237], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1017 OF 15495 ***
    // Wavefunction(s) for diagram number 1017
    // (none)
    // Amplitude(s) for diagram number 1017
    VVV1_0( w_fp[69], w_fp[6], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];

    // *** DIAGRAM 1018 OF 15495 ***
    // Wavefunction(s) for diagram number 1018
    // (none)
    // Amplitude(s) for diagram number 1018
    FFV1_0( w_fp[3], w_fp[236], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1019 OF 15495 ***
    // Wavefunction(s) for diagram number 1019
    VVVV1P0_1( w_fp[66], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[137] );
    VVVV3P0_1( w_fp[66], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[136] );
    VVVV4P0_1( w_fp[66], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[182] );
    // Amplitude(s) for diagram number 1019
    FFV1_0( w_fp[3], w_fp[229], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 1020 OF 15495 ***
    // Wavefunction(s) for diagram number 1020
    // (none)
    // Amplitude(s) for diagram number 1020
    FFV1_0( w_fp[159], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];

    // *** DIAGRAM 1021 OF 15495 ***
    // Wavefunction(s) for diagram number 1021
    // (none)
    // Amplitude(s) for diagram number 1021
    FFV1_0( w_fp[160], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];

    // *** DIAGRAM 1022 OF 15495 ***
    // Wavefunction(s) for diagram number 1022
    FFV1P0_3( w_fp[157], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[75] );
    // Amplitude(s) for diagram number 1022
    VVV1_0( w_fp[68], w_fp[7], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 1023 OF 15495 ***
    // Wavefunction(s) for diagram number 1023
    // (none)
    // Amplitude(s) for diagram number 1023
    FFV1_0( w_fp[160], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1024 OF 15495 ***
    // Wavefunction(s) for diagram number 1024
    // (none)
    // Amplitude(s) for diagram number 1024
    VVV1_0( w_fp[69], w_fp[6], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];

    // *** DIAGRAM 1025 OF 15495 ***
    // Wavefunction(s) for diagram number 1025
    // (none)
    // Amplitude(s) for diagram number 1025
    FFV1_0( w_fp[159], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1026 OF 15495 ***
    // Wavefunction(s) for diagram number 1026
    // (none)
    // Amplitude(s) for diagram number 1026
    FFV1_0( w_fp[157], w_fp[2], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 1027 OF 15495 ***
    // Wavefunction(s) for diagram number 1027
    // (none)
    // Amplitude(s) for diagram number 1027
    FFV1_0( w_fp[166], w_fp[148], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1028 OF 15495 ***
    // Wavefunction(s) for diagram number 1028
    // (none)
    // Amplitude(s) for diagram number 1028
    FFV1_0( w_fp[3], w_fp[148], w_fp[39], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];

    // *** DIAGRAM 1029 OF 15495 ***
    // Wavefunction(s) for diagram number 1029
    // (none)
    // Amplitude(s) for diagram number 1029
    FFV1_0( w_fp[221], w_fp[241], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1030 OF 15495 ***
    // Wavefunction(s) for diagram number 1030
    // (none)
    // Amplitude(s) for diagram number 1030
    FFV1_0( w_fp[221], w_fp[2], w_fp[39], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 1031 OF 15495 ***
    // Wavefunction(s) for diagram number 1031
    // (none)
    // Amplitude(s) for diagram number 1031
    FFV1_0( w_fp[3], w_fp[241], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];

    // *** DIAGRAM 1032 OF 15495 ***
    // Wavefunction(s) for diagram number 1032
    // (none)
    // Amplitude(s) for diagram number 1032
    FFV1_0( w_fp[166], w_fp[2], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 1033 OF 15495 ***
    // Wavefunction(s) for diagram number 1033
    // (none)
    // Amplitude(s) for diagram number 1033
    FFV1_0( w_fp[167], w_fp[148], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[304] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1034 OF 15495 ***
    // Wavefunction(s) for diagram number 1034
    // (none)
    // Amplitude(s) for diagram number 1034
    FFV1_0( w_fp[3], w_fp[148], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];

    // *** DIAGRAM 1035 OF 15495 ***
    // Wavefunction(s) for diagram number 1035
    // (none)
    // Amplitude(s) for diagram number 1035
    FFV1_0( w_fp[221], w_fp[242], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1036 OF 15495 ***
    // Wavefunction(s) for diagram number 1036
    // (none)
    // Amplitude(s) for diagram number 1036
    FFV1_0( w_fp[221], w_fp[2], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 1037 OF 15495 ***
    // Wavefunction(s) for diagram number 1037
    // (none)
    // Amplitude(s) for diagram number 1037
    FFV1_0( w_fp[3], w_fp[242], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1038 OF 15495 ***
    // Wavefunction(s) for diagram number 1038
    // (none)
    // Amplitude(s) for diagram number 1038
    FFV1_0( w_fp[167], w_fp[2], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1039 OF 15495 ***
    // Wavefunction(s) for diagram number 1039
    // (none)
    // Amplitude(s) for diagram number 1039
    FFV1_0( w_fp[3], w_fp[148], w_fp[63], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[148], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[148], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];

    // *** DIAGRAM 1040 OF 15495 ***
    // Wavefunction(s) for diagram number 1040
    // (none)
    // Amplitude(s) for diagram number 1040
    FFV1_0( w_fp[221], w_fp[2], w_fp[63], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    FFV1_0( w_fp[221], w_fp[2], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[221], w_fp[2], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 1041 OF 15495 ***
    // Wavefunction(s) for diagram number 1041
    // (none)
    // Amplitude(s) for diagram number 1041
    FFV1_0( w_fp[221], w_fp[229], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1042 OF 15495 ***
    // Wavefunction(s) for diagram number 1042
    // (none)
    // Amplitude(s) for diagram number 1042
    FFV1_0( w_fp[188], w_fp[229], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1043 OF 15495 ***
    // Wavefunction(s) for diagram number 1043
    VVV1P0_1( w_fp[66], w_fp[84], COUPs[0], 1.0, depCoup, 0., 0., w_fp[65] );
    // Amplitude(s) for diagram number 1043
    FFV1_0( w_fp[3], w_fp[229], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 1044 OF 15495 ***
    // Wavefunction(s) for diagram number 1044
    // (none)
    // Amplitude(s) for diagram number 1044
    FFV1_0( w_fp[157], w_fp[148], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[310] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1045 OF 15495 ***
    // Wavefunction(s) for diagram number 1045
    // (none)
    // Amplitude(s) for diagram number 1045
    FFV1_0( w_fp[157], w_fp[128], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1046 OF 15495 ***
    // Wavefunction(s) for diagram number 1046
    // (none)
    // Amplitude(s) for diagram number 1046
    FFV1_0( w_fp[157], w_fp[2], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 1047 OF 15495 ***
    // Wavefunction(s) for diagram number 1047
    // (none)
    // Amplitude(s) for diagram number 1047
    FFV1_0( w_fp[3], w_fp[128], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 1048 OF 15495 ***
    // Wavefunction(s) for diagram number 1048
    // (none)
    // Amplitude(s) for diagram number 1048
    FFV1_0( w_fp[188], w_fp[2], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[1] += amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];

    // *** DIAGRAM 1049 OF 15495 ***
    // Wavefunction(s) for diagram number 1049
    // (none)
    // Amplitude(s) for diagram number 1049
    FFV1_0( w_fp[3], w_fp[148], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];

    // *** DIAGRAM 1050 OF 15495 ***
    // Wavefunction(s) for diagram number 1050
    // (none)
    // Amplitude(s) for diagram number 1050
    FFV1_0( w_fp[221], w_fp[2], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 1051 OF 15495 ***
    // Wavefunction(s) for diagram number 1051
    // (none)
    // Amplitude(s) for diagram number 1051
    FFV1_0( w_fp[188], w_fp[148], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1052 OF 15495 ***
    // Wavefunction(s) for diagram number 1052
    // (none)
    // Amplitude(s) for diagram number 1052
    FFV1_0( w_fp[221], w_fp[128], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1053 OF 15495 ***
    // Wavefunction(s) for diagram number 1053
    // (none)
    // Amplitude(s) for diagram number 1053
    FFV1_0( w_fp[206], w_fp[238], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];

    // *** DIAGRAM 1054 OF 15495 ***
    // Wavefunction(s) for diagram number 1054
    // (none)
    // Amplitude(s) for diagram number 1054
    FFV1_0( w_fp[206], w_fp[237], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];

    // *** DIAGRAM 1055 OF 15495 ***
    // Wavefunction(s) for diagram number 1055
    // (none)
    // Amplitude(s) for diagram number 1055
    VVV1_0( w_fp[87], w_fp[7], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 1056 OF 15495 ***
    // Wavefunction(s) for diagram number 1056
    // (none)
    // Amplitude(s) for diagram number 1056
    FFV1_0( w_fp[3], w_fp[237], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1057 OF 15495 ***
    // Wavefunction(s) for diagram number 1057
    // (none)
    // Amplitude(s) for diagram number 1057
    VVV1_0( w_fp[88], w_fp[5], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[123] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 1058 OF 15495 ***
    // Wavefunction(s) for diagram number 1058
    // (none)
    // Amplitude(s) for diagram number 1058
    FFV1_0( w_fp[3], w_fp[238], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1059 OF 15495 ***
    // Wavefunction(s) for diagram number 1059
    VVVV1P0_1( w_fp[86], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[85] );
    VVVV3P0_1( w_fp[86], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[70] );
    VVVV4P0_1( w_fp[86], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[64] );
    // Amplitude(s) for diagram number 1059
    FFV1_0( w_fp[3], w_fp[229], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];

    // *** DIAGRAM 1060 OF 15495 ***
    // Wavefunction(s) for diagram number 1060
    // (none)
    // Amplitude(s) for diagram number 1060
    FFV1_0( w_fp[162], w_fp[118], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 1061 OF 15495 ***
    // Wavefunction(s) for diagram number 1061
    // (none)
    // Amplitude(s) for diagram number 1061
    FFV1_0( w_fp[160], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];

    // *** DIAGRAM 1062 OF 15495 ***
    // Wavefunction(s) for diagram number 1062
    // (none)
    // Amplitude(s) for diagram number 1062
    VVV1_0( w_fp[87], w_fp[7], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 1063 OF 15495 ***
    // Wavefunction(s) for diagram number 1063
    // (none)
    // Amplitude(s) for diagram number 1063
    FFV1_0( w_fp[160], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1064 OF 15495 ***
    // Wavefunction(s) for diagram number 1064
    // (none)
    // Amplitude(s) for diagram number 1064
    VVV1_0( w_fp[88], w_fp[5], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 1065 OF 15495 ***
    // Wavefunction(s) for diagram number 1065
    // (none)
    // Amplitude(s) for diagram number 1065
    FFV1_0( w_fp[162], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1066 OF 15495 ***
    // Wavefunction(s) for diagram number 1066
    // (none)
    // Amplitude(s) for diagram number 1066
    FFV1_0( w_fp[157], w_fp[2], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 1067 OF 15495 ***
    // Wavefunction(s) for diagram number 1067
    // (none)
    // Amplitude(s) for diagram number 1067
    FFV1_0( w_fp[165], w_fp[118], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1068 OF 15495 ***
    // Wavefunction(s) for diagram number 1068
    // (none)
    // Amplitude(s) for diagram number 1068
    FFV1_0( w_fp[3], w_fp[118], w_fp[29], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];

    // *** DIAGRAM 1069 OF 15495 ***
    // Wavefunction(s) for diagram number 1069
    // (none)
    // Amplitude(s) for diagram number 1069
    FFV1_0( w_fp[206], w_fp[240], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1070 OF 15495 ***
    // Wavefunction(s) for diagram number 1070
    // (none)
    // Amplitude(s) for diagram number 1070
    FFV1_0( w_fp[206], w_fp[2], w_fp[29], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 1071 OF 15495 ***
    // Wavefunction(s) for diagram number 1071
    // (none)
    // Amplitude(s) for diagram number 1071
    FFV1_0( w_fp[3], w_fp[240], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];

    // *** DIAGRAM 1072 OF 15495 ***
    // Wavefunction(s) for diagram number 1072
    // (none)
    // Amplitude(s) for diagram number 1072
    FFV1_0( w_fp[165], w_fp[2], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 1073 OF 15495 ***
    // Wavefunction(s) for diagram number 1073
    // (none)
    // Amplitude(s) for diagram number 1073
    FFV1_0( w_fp[167], w_fp[118], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1074 OF 15495 ***
    // Wavefunction(s) for diagram number 1074
    // (none)
    // Amplitude(s) for diagram number 1074
    FFV1_0( w_fp[3], w_fp[118], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];

    // *** DIAGRAM 1075 OF 15495 ***
    // Wavefunction(s) for diagram number 1075
    // (none)
    // Amplitude(s) for diagram number 1075
    FFV1_0( w_fp[206], w_fp[242], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1076 OF 15495 ***
    // Wavefunction(s) for diagram number 1076
    // (none)
    // Amplitude(s) for diagram number 1076
    FFV1_0( w_fp[206], w_fp[2], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

    // *** DIAGRAM 1077 OF 15495 ***
    // Wavefunction(s) for diagram number 1077
    // (none)
    // Amplitude(s) for diagram number 1077
    FFV1_0( w_fp[3], w_fp[242], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1078 OF 15495 ***
    // Wavefunction(s) for diagram number 1078
    // (none)
    // Amplitude(s) for diagram number 1078
    FFV1_0( w_fp[167], w_fp[2], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 1079 OF 15495 ***
    // Wavefunction(s) for diagram number 1079
    // (none)
    // Amplitude(s) for diagram number 1079
    FFV1_0( w_fp[3], w_fp[118], w_fp[60], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[118], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[118], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 1080 OF 15495 ***
    // Wavefunction(s) for diagram number 1080
    // (none)
    // Amplitude(s) for diagram number 1080
    FFV1_0( w_fp[206], w_fp[2], w_fp[60], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    FFV1_0( w_fp[206], w_fp[2], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[206], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];

    // *** DIAGRAM 1081 OF 15495 ***
    // Wavefunction(s) for diagram number 1081
    // (none)
    // Amplitude(s) for diagram number 1081
    FFV1_0( w_fp[206], w_fp[229], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1082 OF 15495 ***
    // Wavefunction(s) for diagram number 1082
    // (none)
    // Amplitude(s) for diagram number 1082
    FFV1_0( w_fp[186], w_fp[229], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1083 OF 15495 ***
    // Wavefunction(s) for diagram number 1083
    VVV1P0_1( w_fp[86], w_fp[100], COUPs[0], 1.0, depCoup, 0., 0., w_fp[62] );
    // Amplitude(s) for diagram number 1083
    FFV1_0( w_fp[3], w_fp[229], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];

    // *** DIAGRAM 1084 OF 15495 ***
    // Wavefunction(s) for diagram number 1084
    // (none)
    // Amplitude(s) for diagram number 1084
    FFV1_0( w_fp[157], w_fp[118], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[334] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[544] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1085 OF 15495 ***
    // Wavefunction(s) for diagram number 1085
    // (none)
    // Amplitude(s) for diagram number 1085
    FFV1_0( w_fp[157], w_fp[122], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1086 OF 15495 ***
    // Wavefunction(s) for diagram number 1086
    // (none)
    // Amplitude(s) for diagram number 1086
    FFV1_0( w_fp[157], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 1087 OF 15495 ***
    // Wavefunction(s) for diagram number 1087
    // (none)
    // Amplitude(s) for diagram number 1087
    FFV1_0( w_fp[3], w_fp[122], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 1088 OF 15495 ***
    // Wavefunction(s) for diagram number 1088
    // (none)
    // Amplitude(s) for diagram number 1088
    FFV1_0( w_fp[186], w_fp[2], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];

    // *** DIAGRAM 1089 OF 15495 ***
    // Wavefunction(s) for diagram number 1089
    // (none)
    // Amplitude(s) for diagram number 1089
    FFV1_0( w_fp[3], w_fp[118], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];

    // *** DIAGRAM 1090 OF 15495 ***
    // Wavefunction(s) for diagram number 1090
    // (none)
    // Amplitude(s) for diagram number 1090
    FFV1_0( w_fp[206], w_fp[2], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 1091 OF 15495 ***
    // Wavefunction(s) for diagram number 1091
    // (none)
    // Amplitude(s) for diagram number 1091
    FFV1_0( w_fp[186], w_fp[118], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1092 OF 15495 ***
    // Wavefunction(s) for diagram number 1092
    // (none)
    // Amplitude(s) for diagram number 1092
    FFV1_0( w_fp[206], w_fp[122], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1093 OF 15495 ***
    // Wavefunction(s) for diagram number 1093
    // (none)
    // Amplitude(s) for diagram number 1093
    FFV1_0( w_fp[208], w_fp[238], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];

    // *** DIAGRAM 1094 OF 15495 ***
    // Wavefunction(s) for diagram number 1094
    // (none)
    // Amplitude(s) for diagram number 1094
    FFV1_0( w_fp[208], w_fp[236], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];

    // *** DIAGRAM 1095 OF 15495 ***
    // Wavefunction(s) for diagram number 1095
    // (none)
    // Amplitude(s) for diagram number 1095
    VVV1_0( w_fp[103], w_fp[6], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];

    // *** DIAGRAM 1096 OF 15495 ***
    // Wavefunction(s) for diagram number 1096
    // (none)
    // Amplitude(s) for diagram number 1096
    FFV1_0( w_fp[3], w_fp[236], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1097 OF 15495 ***
    // Wavefunction(s) for diagram number 1097
    // (none)
    // Amplitude(s) for diagram number 1097
    VVV1_0( w_fp[104], w_fp[5], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[125] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];

    // *** DIAGRAM 1098 OF 15495 ***
    // Wavefunction(s) for diagram number 1098
    // (none)
    // Amplitude(s) for diagram number 1098
    FFV1_0( w_fp[3], w_fp[238], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1099 OF 15495 ***
    // Wavefunction(s) for diagram number 1099
    VVVV1P0_1( w_fp[102], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[101] );
    VVVV3P0_1( w_fp[102], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[89] );
    VVVV4P0_1( w_fp[102], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[61] );
    // Amplitude(s) for diagram number 1099
    FFV1_0( w_fp[3], w_fp[229], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];

    // *** DIAGRAM 1100 OF 15495 ***
    // Wavefunction(s) for diagram number 1100
    // (none)
    // Amplitude(s) for diagram number 1100
    FFV1_0( w_fp[162], w_fp[98], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 1101 OF 15495 ***
    // Wavefunction(s) for diagram number 1101
    // (none)
    // Amplitude(s) for diagram number 1101
    FFV1_0( w_fp[159], w_fp[98], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];

    // *** DIAGRAM 1102 OF 15495 ***
    // Wavefunction(s) for diagram number 1102
    // (none)
    // Amplitude(s) for diagram number 1102
    VVV1_0( w_fp[103], w_fp[6], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];

    // *** DIAGRAM 1103 OF 15495 ***
    // Wavefunction(s) for diagram number 1103
    // (none)
    // Amplitude(s) for diagram number 1103
    FFV1_0( w_fp[159], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1104 OF 15495 ***
    // Wavefunction(s) for diagram number 1104
    // (none)
    // Amplitude(s) for diagram number 1104
    VVV1_0( w_fp[104], w_fp[5], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[592] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 1105 OF 15495 ***
    // Wavefunction(s) for diagram number 1105
    // (none)
    // Amplitude(s) for diagram number 1105
    FFV1_0( w_fp[162], w_fp[2], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[550] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1106 OF 15495 ***
    // Wavefunction(s) for diagram number 1106
    // (none)
    // Amplitude(s) for diagram number 1106
    FFV1_0( w_fp[157], w_fp[2], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[89], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];

    // *** DIAGRAM 1107 OF 15495 ***
    // Wavefunction(s) for diagram number 1107
    // (none)
    // Amplitude(s) for diagram number 1107
    FFV1_0( w_fp[165], w_fp[98], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1108 OF 15495 ***
    // Wavefunction(s) for diagram number 1108
    // (none)
    // Amplitude(s) for diagram number 1108
    FFV1_0( w_fp[3], w_fp[98], w_fp[28], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 1109 OF 15495 ***
    // Wavefunction(s) for diagram number 1109
    // (none)
    // Amplitude(s) for diagram number 1109
    FFV1_0( w_fp[208], w_fp[240], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1110 OF 15495 ***
    // Wavefunction(s) for diagram number 1110
    // (none)
    // Amplitude(s) for diagram number 1110
    FFV1_0( w_fp[208], w_fp[2], w_fp[28], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];

    // *** DIAGRAM 1111 OF 15495 ***
    // Wavefunction(s) for diagram number 1111
    // (none)
    // Amplitude(s) for diagram number 1111
    FFV1_0( w_fp[3], w_fp[240], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];

    // *** DIAGRAM 1112 OF 15495 ***
    // Wavefunction(s) for diagram number 1112
    // (none)
    // Amplitude(s) for diagram number 1112
    FFV1_0( w_fp[165], w_fp[2], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 1113 OF 15495 ***
    // Wavefunction(s) for diagram number 1113
    // (none)
    // Amplitude(s) for diagram number 1113
    FFV1_0( w_fp[166], w_fp[98], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1114 OF 15495 ***
    // Wavefunction(s) for diagram number 1114
    // (none)
    // Amplitude(s) for diagram number 1114
    FFV1_0( w_fp[3], w_fp[98], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[356] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];

    // *** DIAGRAM 1115 OF 15495 ***
    // Wavefunction(s) for diagram number 1115
    // (none)
    // Amplitude(s) for diagram number 1115
    FFV1_0( w_fp[208], w_fp[241], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1116 OF 15495 ***
    // Wavefunction(s) for diagram number 1116
    // (none)
    // Amplitude(s) for diagram number 1116
    FFV1_0( w_fp[208], w_fp[2], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];

    // *** DIAGRAM 1117 OF 15495 ***
    // Wavefunction(s) for diagram number 1117
    // (none)
    // Amplitude(s) for diagram number 1117
    FFV1_0( w_fp[3], w_fp[241], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[16] += amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];

    // *** DIAGRAM 1118 OF 15495 ***
    // Wavefunction(s) for diagram number 1118
    // (none)
    // Amplitude(s) for diagram number 1118
    FFV1_0( w_fp[166], w_fp[2], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];

    // *** DIAGRAM 1119 OF 15495 ***
    // Wavefunction(s) for diagram number 1119
    // (none)
    // Amplitude(s) for diagram number 1119
    FFV1_0( w_fp[3], w_fp[98], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[98], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[337] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[98], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[356] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 1120 OF 15495 ***
    // Wavefunction(s) for diagram number 1120
    // (none)
    // Amplitude(s) for diagram number 1120
    FFV1_0( w_fp[208], w_fp[2], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    FFV1_0( w_fp[208], w_fp[2], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    FFV1_0( w_fp[208], w_fp[2], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];

    // *** DIAGRAM 1121 OF 15495 ***
    // Wavefunction(s) for diagram number 1121
    // (none)
    // Amplitude(s) for diagram number 1121
    FFV1_0( w_fp[208], w_fp[229], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1122 OF 15495 ***
    // Wavefunction(s) for diagram number 1122
    // (none)
    // Amplitude(s) for diagram number 1122
    FFV1_0( w_fp[184], w_fp[229], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1123 OF 15495 ***
    // Wavefunction(s) for diagram number 1123
    VVV1P0_1( w_fp[102], w_fp[113], COUPs[0], 1.0, depCoup, 0., 0., w_fp[59] );
    // Amplitude(s) for diagram number 1123
    FFV1_0( w_fp[3], w_fp[229], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];

    // *** DIAGRAM 1124 OF 15495 ***
    // Wavefunction(s) for diagram number 1124
    // (none)
    // Amplitude(s) for diagram number 1124
    FFV1_0( w_fp[157], w_fp[98], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1125 OF 15495 ***
    // Wavefunction(s) for diagram number 1125
    // (none)
    // Amplitude(s) for diagram number 1125
    FFV1_0( w_fp[157], w_fp[244], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1126 OF 15495 ***
    // Wavefunction(s) for diagram number 1126
    // (none)
    // Amplitude(s) for diagram number 1126
    FFV1_0( w_fp[157], w_fp[2], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 1127 OF 15495 ***
    // Wavefunction(s) for diagram number 1127
    // (none)
    // Amplitude(s) for diagram number 1127
    FFV1_0( w_fp[3], w_fp[244], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];

    // *** DIAGRAM 1128 OF 15495 ***
    // Wavefunction(s) for diagram number 1128
    // (none)
    // Amplitude(s) for diagram number 1128
    FFV1_0( w_fp[184], w_fp[2], w_fp[105], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 1129 OF 15495 ***
    // Wavefunction(s) for diagram number 1129
    // (none)
    // Amplitude(s) for diagram number 1129
    FFV1_0( w_fp[3], w_fp[98], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 1130 OF 15495 ***
    // Wavefunction(s) for diagram number 1130
    // (none)
    // Amplitude(s) for diagram number 1130
    FFV1_0( w_fp[208], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];

    // *** DIAGRAM 1131 OF 15495 ***
    // Wavefunction(s) for diagram number 1131
    // (none)
    // Amplitude(s) for diagram number 1131
    FFV1_0( w_fp[184], w_fp[98], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[343] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1132 OF 15495 ***
    // Wavefunction(s) for diagram number 1132
    // (none)
    // Amplitude(s) for diagram number 1132
    FFV1_0( w_fp[208], w_fp[244], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1133 OF 15495 ***
    // Wavefunction(s) for diagram number 1133
    // (none)
    // Amplitude(s) for diagram number 1133
    FFV1_0( w_fp[184], w_fp[155], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];

    // *** DIAGRAM 1134 OF 15495 ***
    // Wavefunction(s) for diagram number 1134
    // (none)
    // Amplitude(s) for diagram number 1134
    FFV1_0( w_fp[184], w_fp[237], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];

    // *** DIAGRAM 1135 OF 15495 ***
    // Wavefunction(s) for diagram number 1135
    // (none)
    // Amplitude(s) for diagram number 1135
    VVV1_0( w_fp[115], w_fp[7], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 1136 OF 15495 ***
    // Wavefunction(s) for diagram number 1136
    // (none)
    // Amplitude(s) for diagram number 1136
    FFV1_0( w_fp[3], w_fp[237], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1137 OF 15495 ***
    // Wavefunction(s) for diagram number 1137
    // (none)
    // Amplitude(s) for diagram number 1137
    VVV1_0( w_fp[4], w_fp[116], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 1138 OF 15495 ***
    // Wavefunction(s) for diagram number 1138
    // (none)
    // Amplitude(s) for diagram number 1138
    FFV1_0( w_fp[3], w_fp[155], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1139 OF 15495 ***
    // Wavefunction(s) for diagram number 1139
    VVVV1P0_1( w_fp[4], w_fp[113], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[237] );
    VVVV3P0_1( w_fp[4], w_fp[113], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[8] );
    VVVV4P0_1( w_fp[4], w_fp[113], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[114] );
    // Amplitude(s) for diagram number 1139
    FFV1_0( w_fp[3], w_fp[229], w_fp[237], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 1140 OF 15495 ***
    // Wavefunction(s) for diagram number 1140
    // (none)
    // Amplitude(s) for diagram number 1140
    FFV1_0( w_fp[192], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 1141 OF 15495 ***
    // Wavefunction(s) for diagram number 1141
    // (none)
    // Amplitude(s) for diagram number 1141
    FFV1_0( w_fp[160], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];

    // *** DIAGRAM 1142 OF 15495 ***
    // Wavefunction(s) for diagram number 1142
    // (none)
    // Amplitude(s) for diagram number 1142
    VVV1_0( w_fp[115], w_fp[7], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 1143 OF 15495 ***
    // Wavefunction(s) for diagram number 1143
    // (none)
    // Amplitude(s) for diagram number 1143
    FFV1_0( w_fp[160], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[328] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1144 OF 15495 ***
    // Wavefunction(s) for diagram number 1144
    // (none)
    // Amplitude(s) for diagram number 1144
    VVV1_0( w_fp[4], w_fp[116], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 1145 OF 15495 ***
    // Wavefunction(s) for diagram number 1145
    // (none)
    // Amplitude(s) for diagram number 1145
    FFV1_0( w_fp[192], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1146 OF 15495 ***
    // Wavefunction(s) for diagram number 1146
    // (none)
    // Amplitude(s) for diagram number 1146
    FFV1_0( w_fp[157], w_fp[2], w_fp[237], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 1147 OF 15495 ***
    // Wavefunction(s) for diagram number 1147
    // (none)
    // Amplitude(s) for diagram number 1147
    FFV1_0( w_fp[195], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1148 OF 15495 ***
    // Wavefunction(s) for diagram number 1148
    // (none)
    // Amplitude(s) for diagram number 1148
    FFV1_0( w_fp[3], w_fp[244], w_fp[13], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];

    // *** DIAGRAM 1149 OF 15495 ***
    // Wavefunction(s) for diagram number 1149
    // (none)
    // Amplitude(s) for diagram number 1149
    FFV1_0( w_fp[184], w_fp[67], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1150 OF 15495 ***
    // Wavefunction(s) for diagram number 1150
    // (none)
    // Amplitude(s) for diagram number 1150
    FFV1_0( w_fp[184], w_fp[2], w_fp[13], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 1151 OF 15495 ***
    // Wavefunction(s) for diagram number 1151
    // (none)
    // Amplitude(s) for diagram number 1151
    FFV1_0( w_fp[3], w_fp[67], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[264] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];

    // *** DIAGRAM 1152 OF 15495 ***
    // Wavefunction(s) for diagram number 1152
    // (none)
    // Amplitude(s) for diagram number 1152
    FFV1_0( w_fp[195], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 1153 OF 15495 ***
    // Wavefunction(s) for diagram number 1153
    // (none)
    // Amplitude(s) for diagram number 1153
    FFV1_0( w_fp[167], w_fp[244], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[568] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[569] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1154 OF 15495 ***
    // Wavefunction(s) for diagram number 1154
    // (none)
    // Amplitude(s) for diagram number 1154
    FFV1_0( w_fp[3], w_fp[244], w_fp[47], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];

    // *** DIAGRAM 1155 OF 15495 ***
    // Wavefunction(s) for diagram number 1155
    // (none)
    // Amplitude(s) for diagram number 1155
    FFV1_0( w_fp[184], w_fp[242], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1156 OF 15495 ***
    // Wavefunction(s) for diagram number 1156
    // (none)
    // Amplitude(s) for diagram number 1156
    FFV1_0( w_fp[184], w_fp[2], w_fp[47], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[343] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];

    // *** DIAGRAM 1157 OF 15495 ***
    // Wavefunction(s) for diagram number 1157
    // (none)
    // Amplitude(s) for diagram number 1157
    FFV1_0( w_fp[3], w_fp[242], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[19] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 1158 OF 15495 ***
    // Wavefunction(s) for diagram number 1158
    // (none)
    // Amplitude(s) for diagram number 1158
    FFV1_0( w_fp[167], w_fp[2], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[304] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[328] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];

    // *** DIAGRAM 1159 OF 15495 ***
    // Wavefunction(s) for diagram number 1159
    // (none)
    // Amplitude(s) for diagram number 1159
    FFV1_0( w_fp[3], w_fp[244], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[569] += amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[244], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[3], w_fp[244], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[552] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];

    // *** DIAGRAM 1160 OF 15495 ***
    // Wavefunction(s) for diagram number 1160
    // (none)
    // Amplitude(s) for diagram number 1160
    FFV1_0( w_fp[184], w_fp[2], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[18] += amp_sv[0];
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[343] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    FFV1_0( w_fp[184], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[184], w_fp[2], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 1161 OF 15495 ***
    // Wavefunction(s) for diagram number 1161
    // (none)
    // Amplitude(s) for diagram number 1161
    FFV1_0( w_fp[186], w_fp[155], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];

    // *** DIAGRAM 1162 OF 15495 ***
    // Wavefunction(s) for diagram number 1162
    // (none)
    // Amplitude(s) for diagram number 1162
    FFV1_0( w_fp[186], w_fp[236], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];

    // *** DIAGRAM 1163 OF 15495 ***
    // Wavefunction(s) for diagram number 1163
    // (none)
    // Amplitude(s) for diagram number 1163
    VVV1_0( w_fp[124], w_fp[6], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];

    // *** DIAGRAM 1164 OF 15495 ***
    // Wavefunction(s) for diagram number 1164
    // (none)
    // Amplitude(s) for diagram number 1164
    FFV1_0( w_fp[3], w_fp[236], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1165 OF 15495 ***
    // Wavefunction(s) for diagram number 1165
    // (none)
    // Amplitude(s) for diagram number 1165
    VVV1_0( w_fp[4], w_fp[125], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[121] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];

    // *** DIAGRAM 1166 OF 15495 ***
    // Wavefunction(s) for diagram number 1166
    // (none)
    // Amplitude(s) for diagram number 1166
    FFV1_0( w_fp[3], w_fp[155], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[2] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1167 OF 15495 ***
    // Wavefunction(s) for diagram number 1167
    VVVV1P0_1( w_fp[4], w_fp[100], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[236] );
    VVVV3P0_1( w_fp[4], w_fp[100], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[56] );
    VVVV4P0_1( w_fp[4], w_fp[100], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[55] );
    // Amplitude(s) for diagram number 1167
    FFV1_0( w_fp[3], w_fp[229], w_fp[236], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[121] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];

    // *** DIAGRAM 1168 OF 15495 ***
    // Wavefunction(s) for diagram number 1168
    // (none)
    // Amplitude(s) for diagram number 1168
    FFV1_0( w_fp[192], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 1169 OF 15495 ***
    // Wavefunction(s) for diagram number 1169
    // (none)
    // Amplitude(s) for diagram number 1169
    FFV1_0( w_fp[159], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];

    // *** DIAGRAM 1170 OF 15495 ***
    // Wavefunction(s) for diagram number 1170
    // (none)
    // Amplitude(s) for diagram number 1170
    VVV1_0( w_fp[124], w_fp[6], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[550] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];

    // *** DIAGRAM 1171 OF 15495 ***
    // Wavefunction(s) for diagram number 1171
    // (none)
    // Amplitude(s) for diagram number 1171
    FFV1_0( w_fp[159], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[352] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1172 OF 15495 ***
    // Wavefunction(s) for diagram number 1172
    // (none)
    // Amplitude(s) for diagram number 1172
    VVV1_0( w_fp[4], w_fp[125], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 1173 OF 15495 ***
    // Wavefunction(s) for diagram number 1173
    // (none)
    // Amplitude(s) for diagram number 1173
    FFV1_0( w_fp[192], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1174 OF 15495 ***
    // Wavefunction(s) for diagram number 1174
    // (none)
    // Amplitude(s) for diagram number 1174
    FFV1_0( w_fp[157], w_fp[2], w_fp[236], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] -= amp_sv[0];
    jamp_sv[311] += amp_sv[0];
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[352] += amp_sv[0];
    jamp_sv[353] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[56], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[328] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    FFV1_0( w_fp[157], w_fp[2], w_fp[55], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[550] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 1175 OF 15495 ***
    // Wavefunction(s) for diagram number 1175
    // (none)
    // Amplitude(s) for diagram number 1175
    FFV1_0( w_fp[195], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1176 OF 15495 ***
    // Wavefunction(s) for diagram number 1176
    // (none)
    // Amplitude(s) for diagram number 1176
    FFV1_0( w_fp[3], w_fp[122], w_fp[12], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];

    // *** DIAGRAM 1177 OF 15495 ***
    // Wavefunction(s) for diagram number 1177
    // (none)
    // Amplitude(s) for diagram number 1177
    FFV1_0( w_fp[186], w_fp[67], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1178 OF 15495 ***
    // Wavefunction(s) for diagram number 1178
    // (none)
    // Amplitude(s) for diagram number 1178
    FFV1_0( w_fp[186], w_fp[2], w_fp[12], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];

    // *** DIAGRAM 1179 OF 15495 ***
    // Wavefunction(s) for diagram number 1179
    // (none)
    // Amplitude(s) for diagram number 1179
    FFV1_0( w_fp[3], w_fp[67], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];

    // *** DIAGRAM 1180 OF 15495 ***
    // Wavefunction(s) for diagram number 1180
    // (none)
    // Amplitude(s) for diagram number 1180
    FFV1_0( w_fp[195], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 1181 OF 15495 ***
    // Wavefunction(s) for diagram number 1181
    // (none)
    // Amplitude(s) for diagram number 1181
    FFV1_0( w_fp[166], w_fp[122], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1182 OF 15495 ***
    // Wavefunction(s) for diagram number 1182
    // (none)
    // Amplitude(s) for diagram number 1182
    FFV1_0( w_fp[3], w_fp[122], w_fp[40], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 1183 OF 15495 ***
    // Wavefunction(s) for diagram number 1183
    // (none)
    // Amplitude(s) for diagram number 1183
    FFV1_0( w_fp[186], w_fp[241], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1184 OF 15495 ***
    // Wavefunction(s) for diagram number 1184
    // (none)
    // Amplitude(s) for diagram number 1184
    FFV1_0( w_fp[186], w_fp[2], w_fp[40], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];

    // *** DIAGRAM 1185 OF 15495 ***
    // Wavefunction(s) for diagram number 1185
    // (none)
    // Amplitude(s) for diagram number 1185
    FFV1_0( w_fp[3], w_fp[241], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[13] += amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];

    // *** DIAGRAM 1186 OF 15495 ***
    // Wavefunction(s) for diagram number 1186
    // (none)
    // Amplitude(s) for diagram number 1186
    FFV1_0( w_fp[166], w_fp[2], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[352] -= amp_sv[0];
    jamp_sv[353] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];

    // *** DIAGRAM 1187 OF 15495 ***
    // Wavefunction(s) for diagram number 1187
    // (none)
    // Amplitude(s) for diagram number 1187
    FFV1_0( w_fp[3], w_fp[122], w_fp[51], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[122], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[689] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[122], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 1188 OF 15495 ***
    // Wavefunction(s) for diagram number 1188
    // (none)
    // Amplitude(s) for diagram number 1188
    FFV1_0( w_fp[186], w_fp[2], w_fp[51], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[3] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    FFV1_0( w_fp[186], w_fp[2], w_fp[52], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[13] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    FFV1_0( w_fp[186], w_fp[2], w_fp[53], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 1189 OF 15495 ***
    // Wavefunction(s) for diagram number 1189
    // (none)
    // Amplitude(s) for diagram number 1189
    FFV1_0( w_fp[188], w_fp[155], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];

    // *** DIAGRAM 1190 OF 15495 ***
    // Wavefunction(s) for diagram number 1190
    // (none)
    // Amplitude(s) for diagram number 1190
    FFV1_0( w_fp[188], w_fp[238], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];

    // *** DIAGRAM 1191 OF 15495 ***
    // Wavefunction(s) for diagram number 1191
    // (none)
    // Amplitude(s) for diagram number 1191
    VVV1_0( w_fp[130], w_fp[5], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];

    // *** DIAGRAM 1192 OF 15495 ***
    // Wavefunction(s) for diagram number 1192
    // (none)
    // Amplitude(s) for diagram number 1192
    FFV1_0( w_fp[3], w_fp[238], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1193 OF 15495 ***
    // Wavefunction(s) for diagram number 1193
    // (none)
    // Amplitude(s) for diagram number 1193
    VVV1_0( w_fp[4], w_fp[131], w_fp[138], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];

    // *** DIAGRAM 1194 OF 15495 ***
    // Wavefunction(s) for diagram number 1194
    // (none)
    // Amplitude(s) for diagram number 1194
    FFV1_0( w_fp[3], w_fp[155], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1195 OF 15495 ***
    // Wavefunction(s) for diagram number 1195
    VVVV1P0_1( w_fp[4], w_fp[5], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[155] );
    VVVV3P0_1( w_fp[4], w_fp[5], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[138] );
    VVVV4P0_1( w_fp[4], w_fp[5], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[238] );
    // Amplitude(s) for diagram number 1195
    FFV1_0( w_fp[3], w_fp[229], w_fp[155], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[138], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[16] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 1196 OF 15495 ***
    // Wavefunction(s) for diagram number 1196
    // (none)
    // Amplitude(s) for diagram number 1196
    FFV1_0( w_fp[192], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[598] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 1197 OF 15495 ***
    // Wavefunction(s) for diagram number 1197
    // (none)
    // Amplitude(s) for diagram number 1197
    FFV1_0( w_fp[162], w_fp[128], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 1198 OF 15495 ***
    // Wavefunction(s) for diagram number 1198
    // (none)
    // Amplitude(s) for diagram number 1198
    VVV1_0( w_fp[130], w_fp[5], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[358] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    jamp_sv[592] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 1199 OF 15495 ***
    // Wavefunction(s) for diagram number 1199
    // (none)
    // Amplitude(s) for diagram number 1199
    FFV1_0( w_fp[162], w_fp[2], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[334] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[358] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[592] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 1200 OF 15495 ***
    // Wavefunction(s) for diagram number 1200
    // (none)
    // Amplitude(s) for diagram number 1200
    VVV1_0( w_fp[4], w_fp[131], w_fp[75], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[304] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[310] += amp_sv[0];
    jamp_sv[311] -= amp_sv[0];
    jamp_sv[334] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[358] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 8 );
    storeWf( wfs, w_cx, nevt, 55 );
    storeWf( wfs, w_cx, nevt, 56 );
    storeWf( wfs, w_cx, nevt, 59 );
    storeWf( wfs, w_cx, nevt, 61 );
    storeWf( wfs, w_cx, nevt, 62 );
    storeWf( wfs, w_cx, nevt, 64 );
    storeWf( wfs, w_cx, nevt, 65 );
    storeWf( wfs, w_cx, nevt, 70 );
    storeWf( wfs, w_cx, nevt, 75 );
    storeWf( wfs, w_cx, nevt, 85 );
    storeWf( wfs, w_cx, nevt, 89 );
    storeWf( wfs, w_cx, nevt, 101 );
    storeWf( wfs, w_cx, nevt, 114 );
    storeWf( wfs, w_cx, nevt, 136 );
    storeWf( wfs, w_cx, nevt, 137 );
    storeWf( wfs, w_cx, nevt, 138 );
    storeWf( wfs, w_cx, nevt, 155 );
    storeWf( wfs, w_cx, nevt, 182 );
    storeWf( wfs, w_cx, nevt, 236 );
    storeWf( wfs, w_cx, nevt, 237 );
    storeWf( wfs, w_cx, nevt, 238 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
