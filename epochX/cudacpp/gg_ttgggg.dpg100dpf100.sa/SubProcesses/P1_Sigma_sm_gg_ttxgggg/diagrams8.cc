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
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 28 );
    retrieveWf( wfs, w_cx, nevt, 29 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 38 );
    retrieveWf( wfs, w_cx, nevt, 39 );
    retrieveWf( wfs, w_cx, nevt, 44 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 57 );
    retrieveWf( wfs, w_cx, nevt, 58 );
    retrieveWf( wfs, w_cx, nevt, 59 );
    retrieveWf( wfs, w_cx, nevt, 60 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 63 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 67 );
    retrieveWf( wfs, w_cx, nevt, 75 );
    retrieveWf( wfs, w_cx, nevt, 76 );
    retrieveWf( wfs, w_cx, nevt, 77 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 101 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 114 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 119 );
    retrieveWf( wfs, w_cx, nevt, 120 );
    retrieveWf( wfs, w_cx, nevt, 121 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 127 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 138 );
    retrieveWf( wfs, w_cx, nevt, 151 );
    retrieveWf( wfs, w_cx, nevt, 152 );
    retrieveWf( wfs, w_cx, nevt, 153 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 155 );
    retrieveWf( wfs, w_cx, nevt, 157 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 229 );
    retrieveWf( wfs, w_cx, nevt, 235 );
#endif
#endif

    // *** DIAGRAM 701 OF 15495 ***
    // Wavefunction(s) for diagram number 701
    // (none)
    // Amplitude(s) for diagram number 701
    FFV1_0( w_fp[157], w_fp[215], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 702 OF 15495 ***
    // Wavefunction(s) for diagram number 702
    // (none)
    // Amplitude(s) for diagram number 702
    VVV1_0( w_fp[10], w_fp[228], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 703 OF 15495 ***
    // Wavefunction(s) for diagram number 703
    // (none)
    // Amplitude(s) for diagram number 703
    FFV1_0( w_fp[3], w_fp[235], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 704 OF 15495 ***
    // Wavefunction(s) for diagram number 704
    // (none)
    // Amplitude(s) for diagram number 704
    FFV1_0( w_fp[184], w_fp[215], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 705 OF 15495 ***
    // Wavefunction(s) for diagram number 705
    // (none)
    // Amplitude(s) for diagram number 705
    VVV1_0( w_fp[114], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 706 OF 15495 ***
    // Wavefunction(s) for diagram number 706
    // (none)
    // Amplitude(s) for diagram number 706
    FFV1_0( w_fp[3], w_fp[225], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 707 OF 15495 ***
    // Wavefunction(s) for diagram number 707
    // (none)
    // Amplitude(s) for diagram number 707
    VVV1_0( w_fp[8], w_fp[228], w_fp[115], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 708 OF 15495 ***
    // Wavefunction(s) for diagram number 708
    // (none)
    // Amplitude(s) for diagram number 708
    FFV1_0( w_fp[184], w_fp[225], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 709 OF 15495 ***
    // Wavefunction(s) for diagram number 709
    // (none)
    // Amplitude(s) for diagram number 709
    FFV1_0( w_fp[3], w_fp[215], w_fp[77], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[76], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[75], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 710 OF 15495 ***
    // Wavefunction(s) for diagram number 710
    // (none)
    // Amplitude(s) for diagram number 710
    FFV1_0( w_fp[3], w_fp[229], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[229], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 711 OF 15495 ***
    // Wavefunction(s) for diagram number 711
    // (none)
    // Amplitude(s) for diagram number 711
    FFV1_0( w_fp[157], w_fp[215], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[157], w_fp[215], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 712 OF 15495 ***
    // Wavefunction(s) for diagram number 712
    // (none)
    // Amplitude(s) for diagram number 712
    FFV1_0( w_fp[3], w_fp[215], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[665] -= amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[718] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[138], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[665] += amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[689] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[718] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 713 OF 15495 ***
    // Wavefunction(s) for diagram number 713
    FFV1_1( w_fp[2], w_fp[8], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[229] );
    FFV1_1( w_fp[229], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[236] );
    // Amplitude(s) for diagram number 713
    FFV1_0( w_fp[216], w_fp[236], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 714 OF 15495 ***
    // Wavefunction(s) for diagram number 714
    FFV1_1( w_fp[229], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[237] );
    // Amplitude(s) for diagram number 714
    FFV1_0( w_fp[216], w_fp[237], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 715 OF 15495 ***
    // Wavefunction(s) for diagram number 715
    FFV1_1( w_fp[229], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[238] );
    // Amplitude(s) for diagram number 715
    FFV1_0( w_fp[198], w_fp[238], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 716 OF 15495 ***
    // Wavefunction(s) for diagram number 716
    // (none)
    // Amplitude(s) for diagram number 716
    FFV1_0( w_fp[198], w_fp[237], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 717 OF 15495 ***
    // Wavefunction(s) for diagram number 717
    // (none)
    // Amplitude(s) for diagram number 717
    FFV1_0( w_fp[199], w_fp[238], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 718 OF 15495 ***
    // Wavefunction(s) for diagram number 718
    // (none)
    // Amplitude(s) for diagram number 718
    FFV1_0( w_fp[199], w_fp[236], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 719 OF 15495 ***
    // Wavefunction(s) for diagram number 719
    FFV1P0_3( w_fp[196], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[239] );
    // Amplitude(s) for diagram number 719
    VVVV1_0( w_fp[26], w_fp[239], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV3_0( w_fp[26], w_fp[239], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVVV4_0( w_fp[26], w_fp[239], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 720 OF 15495 ***
    // Wavefunction(s) for diagram number 720
    // (none)
    // Amplitude(s) for diagram number 720
    VVV1_0( w_fp[239], w_fp[7], w_fp[28], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 721 OF 15495 ***
    // Wavefunction(s) for diagram number 721
    // (none)
    // Amplitude(s) for diagram number 721
    VVV1_0( w_fp[239], w_fp[6], w_fp[29], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 722 OF 15495 ***
    // Wavefunction(s) for diagram number 722
    FFV1_1( w_fp[2], w_fp[26], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[240] );
    // Amplitude(s) for diagram number 722
    FFV1_0( w_fp[198], w_fp[240], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 723 OF 15495 ***
    // Wavefunction(s) for diagram number 723
    // (none)
    // Amplitude(s) for diagram number 723
    FFV1_0( w_fp[198], w_fp[2], w_fp[29], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 724 OF 15495 ***
    // Wavefunction(s) for diagram number 724
    // (none)
    // Amplitude(s) for diagram number 724
    FFV1_0( w_fp[199], w_fp[240], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];

    // *** DIAGRAM 725 OF 15495 ***
    // Wavefunction(s) for diagram number 725
    // (none)
    // Amplitude(s) for diagram number 725
    FFV1_0( w_fp[199], w_fp[2], w_fp[28], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 726 OF 15495 ***
    // Wavefunction(s) for diagram number 726
    // (none)
    // Amplitude(s) for diagram number 726
    VVVV1_0( w_fp[37], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    VVVV3_0( w_fp[37], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVVV4_0( w_fp[37], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 727 OF 15495 ***
    // Wavefunction(s) for diagram number 727
    // (none)
    // Amplitude(s) for diagram number 727
    VVV1_0( w_fp[239], w_fp[7], w_fp[38], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 728 OF 15495 ***
    // Wavefunction(s) for diagram number 728
    // (none)
    // Amplitude(s) for diagram number 728
    VVV1_0( w_fp[239], w_fp[5], w_fp[39], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 729 OF 15495 ***
    // Wavefunction(s) for diagram number 729
    FFV1_1( w_fp[2], w_fp[37], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[241] );
    // Amplitude(s) for diagram number 729
    FFV1_0( w_fp[216], w_fp[241], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];

    // *** DIAGRAM 730 OF 15495 ***
    // Wavefunction(s) for diagram number 730
    // (none)
    // Amplitude(s) for diagram number 730
    FFV1_0( w_fp[216], w_fp[2], w_fp[39], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 731 OF 15495 ***
    // Wavefunction(s) for diagram number 731
    // (none)
    // Amplitude(s) for diagram number 731
    FFV1_0( w_fp[199], w_fp[241], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];

    // *** DIAGRAM 732 OF 15495 ***
    // Wavefunction(s) for diagram number 732
    // (none)
    // Amplitude(s) for diagram number 732
    FFV1_0( w_fp[199], w_fp[2], w_fp[38], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 733 OF 15495 ***
    // Wavefunction(s) for diagram number 733
    // (none)
    // Amplitude(s) for diagram number 733
    VVVV1_0( w_fp[44], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    VVVV3_0( w_fp[44], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    VVVV4_0( w_fp[44], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 734 OF 15495 ***
    // Wavefunction(s) for diagram number 734
    // (none)
    // Amplitude(s) for diagram number 734
    VVV1_0( w_fp[239], w_fp[6], w_fp[45], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

    // *** DIAGRAM 735 OF 15495 ***
    // Wavefunction(s) for diagram number 735
    // (none)
    // Amplitude(s) for diagram number 735
    VVV1_0( w_fp[239], w_fp[5], w_fp[46], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 736 OF 15495 ***
    // Wavefunction(s) for diagram number 736
    FFV1_1( w_fp[2], w_fp[44], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[242] );
    // Amplitude(s) for diagram number 736
    FFV1_0( w_fp[216], w_fp[242], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 737 OF 15495 ***
    // Wavefunction(s) for diagram number 737
    // (none)
    // Amplitude(s) for diagram number 737
    FFV1_0( w_fp[216], w_fp[2], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 738 OF 15495 ***
    // Wavefunction(s) for diagram number 738
    // (none)
    // Amplitude(s) for diagram number 738
    FFV1_0( w_fp[198], w_fp[242], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

    // *** DIAGRAM 739 OF 15495 ***
    // Wavefunction(s) for diagram number 739
    // (none)
    // Amplitude(s) for diagram number 739
    FFV1_0( w_fp[198], w_fp[2], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 740 OF 15495 ***
    // Wavefunction(s) for diagram number 740
    // (none)
    // Amplitude(s) for diagram number 740
    VVV1_0( w_fp[57], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    VVV1_0( w_fp[58], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    VVV1_0( w_fp[59], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 741 OF 15495 ***
    // Wavefunction(s) for diagram number 741
    // (none)
    // Amplitude(s) for diagram number 741
    FFV1_0( w_fp[199], w_fp[2], w_fp[57], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[58], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[59], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 742 OF 15495 ***
    // Wavefunction(s) for diagram number 742
    // (none)
    // Amplitude(s) for diagram number 742
    VVV1_0( w_fp[60], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    VVV1_0( w_fp[61], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    VVV1_0( w_fp[62], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];

    // *** DIAGRAM 743 OF 15495 ***
    // Wavefunction(s) for diagram number 743
    // (none)
    // Amplitude(s) for diagram number 743
    FFV1_0( w_fp[198], w_fp[2], w_fp[60], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[62], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 744 OF 15495 ***
    // Wavefunction(s) for diagram number 744
    // (none)
    // Amplitude(s) for diagram number 744
    VVV1_0( w_fp[63], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    VVV1_0( w_fp[64], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    VVV1_0( w_fp[65], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 745 OF 15495 ***
    // Wavefunction(s) for diagram number 745
    // (none)
    // Amplitude(s) for diagram number 745
    FFV1_0( w_fp[216], w_fp[2], w_fp[63], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 746 OF 15495 ***
    // Wavefunction(s) for diagram number 746
    FFV1_2( w_fp[196], w_fp[113], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[243] );
    // Amplitude(s) for diagram number 746
    FFV1_0( w_fp[243], w_fp[229], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 747 OF 15495 ***
    // Wavefunction(s) for diagram number 747
    // (none)
    // Amplitude(s) for diagram number 747
    FFV1_0( w_fp[199], w_fp[229], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];

    // *** DIAGRAM 748 OF 15495 ***
    // Wavefunction(s) for diagram number 748
    // (none)
    // Amplitude(s) for diagram number 748
    FFV1_0( w_fp[196], w_fp[229], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 749 OF 15495 ***
    // Wavefunction(s) for diagram number 749
    FFV1_1( w_fp[2], w_fp[113], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[244] );
    // Amplitude(s) for diagram number 749
    FFV1_0( w_fp[200], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];

    // *** DIAGRAM 750 OF 15495 ***
    // Wavefunction(s) for diagram number 750
    // (none)
    // Amplitude(s) for diagram number 750
    FFV1_0( w_fp[200], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 751 OF 15495 ***
    // Wavefunction(s) for diagram number 751
    // (none)
    // Amplitude(s) for diagram number 751
    VVV1_0( w_fp[114], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 752 OF 15495 ***
    // Wavefunction(s) for diagram number 752
    // (none)
    // Amplitude(s) for diagram number 752
    FFV1_0( w_fp[199], w_fp[2], w_fp[114], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 753 OF 15495 ***
    // Wavefunction(s) for diagram number 753
    // (none)
    // Amplitude(s) for diagram number 753
    VVV1_0( w_fp[44], w_fp[239], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 754 OF 15495 ***
    // Wavefunction(s) for diagram number 754
    // (none)
    // Amplitude(s) for diagram number 754
    FFV1_0( w_fp[196], w_fp[244], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 755 OF 15495 ***
    // Wavefunction(s) for diagram number 755
    // (none)
    // Amplitude(s) for diagram number 755
    FFV1_0( w_fp[243], w_fp[2], w_fp[44], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 756 OF 15495 ***
    // Wavefunction(s) for diagram number 756
    // (none)
    // Amplitude(s) for diagram number 756
    VVV1_0( w_fp[8], w_fp[239], w_fp[116], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 757 OF 15495 ***
    // Wavefunction(s) for diagram number 757
    // (none)
    // Amplitude(s) for diagram number 757
    FFV1_0( w_fp[199], w_fp[244], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[433] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];

    // *** DIAGRAM 758 OF 15495 ***
    // Wavefunction(s) for diagram number 758
    // (none)
    // Amplitude(s) for diagram number 758
    FFV1_0( w_fp[196], w_fp[2], w_fp[121], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[122], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[21] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[123], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 759 OF 15495 ***
    // Wavefunction(s) for diagram number 759
    FFV1_2( w_fp[196], w_fp[100], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[123] );
    // Amplitude(s) for diagram number 759
    FFV1_0( w_fp[123], w_fp[229], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];

    // *** DIAGRAM 760 OF 15495 ***
    // Wavefunction(s) for diagram number 760
    // (none)
    // Amplitude(s) for diagram number 760
    FFV1_0( w_fp[198], w_fp[229], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];

    // *** DIAGRAM 761 OF 15495 ***
    // Wavefunction(s) for diagram number 761
    // (none)
    // Amplitude(s) for diagram number 761
    FFV1_0( w_fp[196], w_fp[229], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 762 OF 15495 ***
    // Wavefunction(s) for diagram number 762
    FFV1_1( w_fp[2], w_fp[100], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[122] );
    // Amplitude(s) for diagram number 762
    FFV1_0( w_fp[200], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 763 OF 15495 ***
    // Wavefunction(s) for diagram number 763
    // (none)
    // Amplitude(s) for diagram number 763
    FFV1_0( w_fp[200], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 764 OF 15495 ***
    // Wavefunction(s) for diagram number 764
    // (none)
    // Amplitude(s) for diagram number 764
    VVV1_0( w_fp[101], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 765 OF 15495 ***
    // Wavefunction(s) for diagram number 765
    // (none)
    // Amplitude(s) for diagram number 765
    FFV1_0( w_fp[198], w_fp[2], w_fp[101], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 766 OF 15495 ***
    // Wavefunction(s) for diagram number 766
    // (none)
    // Amplitude(s) for diagram number 766
    VVV1_0( w_fp[37], w_fp[239], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];

    // *** DIAGRAM 767 OF 15495 ***
    // Wavefunction(s) for diagram number 767
    // (none)
    // Amplitude(s) for diagram number 767
    FFV1_0( w_fp[196], w_fp[122], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 768 OF 15495 ***
    // Wavefunction(s) for diagram number 768
    // (none)
    // Amplitude(s) for diagram number 768
    FFV1_0( w_fp[123], w_fp[2], w_fp[37], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 769 OF 15495 ***
    // Wavefunction(s) for diagram number 769
    // (none)
    // Amplitude(s) for diagram number 769
    VVV1_0( w_fp[8], w_fp[239], w_fp[125], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];

    // *** DIAGRAM 770 OF 15495 ***
    // Wavefunction(s) for diagram number 770
    // (none)
    // Amplitude(s) for diagram number 770
    FFV1_0( w_fp[198], w_fp[122], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[457] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];

    // *** DIAGRAM 771 OF 15495 ***
    // Wavefunction(s) for diagram number 771
    // (none)
    // Amplitude(s) for diagram number 771
    FFV1_0( w_fp[196], w_fp[2], w_fp[127], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[128], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[129], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 772 OF 15495 ***
    // Wavefunction(s) for diagram number 772
    // (none)
    // Amplitude(s) for diagram number 772
    FFV1_0( w_fp[216], w_fp[229], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 773 OF 15495 ***
    // Wavefunction(s) for diagram number 773
    FFV1_2( w_fp[196], w_fp[84], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[129] );
    // Amplitude(s) for diagram number 773
    FFV1_0( w_fp[129], w_fp[229], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];

    // *** DIAGRAM 774 OF 15495 ***
    // Wavefunction(s) for diagram number 774
    // (none)
    // Amplitude(s) for diagram number 774
    FFV1_0( w_fp[196], w_fp[229], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 775 OF 15495 ***
    // Wavefunction(s) for diagram number 775
    FFV1_1( w_fp[2], w_fp[84], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[128] );
    // Amplitude(s) for diagram number 775
    FFV1_0( w_fp[200], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 776 OF 15495 ***
    // Wavefunction(s) for diagram number 776
    // (none)
    // Amplitude(s) for diagram number 776
    FFV1_0( w_fp[200], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 777 OF 15495 ***
    // Wavefunction(s) for diagram number 777
    // (none)
    // Amplitude(s) for diagram number 777
    VVV1_0( w_fp[26], w_fp[239], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 778 OF 15495 ***
    // Wavefunction(s) for diagram number 778
    // (none)
    // Amplitude(s) for diagram number 778
    FFV1_0( w_fp[196], w_fp[128], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 779 OF 15495 ***
    // Wavefunction(s) for diagram number 779
    // (none)
    // Amplitude(s) for diagram number 779
    FFV1_0( w_fp[129], w_fp[2], w_fp[26], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 780 OF 15495 ***
    // Wavefunction(s) for diagram number 780
    // (none)
    // Amplitude(s) for diagram number 780
    VVV1_0( w_fp[85], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];

    // *** DIAGRAM 781 OF 15495 ***
    // Wavefunction(s) for diagram number 781
    // (none)
    // Amplitude(s) for diagram number 781
    FFV1_0( w_fp[216], w_fp[2], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 782 OF 15495 ***
    // Wavefunction(s) for diagram number 782
    // (none)
    // Amplitude(s) for diagram number 782
    VVV1_0( w_fp[8], w_fp[239], w_fp[131], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[17] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 783 OF 15495 ***
    // Wavefunction(s) for diagram number 783
    // (none)
    // Amplitude(s) for diagram number 783
    FFV1_0( w_fp[216], w_fp[128], w_fp[8], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];

    // *** DIAGRAM 784 OF 15495 ***
    // Wavefunction(s) for diagram number 784
    // (none)
    // Amplitude(s) for diagram number 784
    FFV1_0( w_fp[196], w_fp[2], w_fp[120], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[119], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[118], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 785 OF 15495 ***
    // Wavefunction(s) for diagram number 785
    // (none)
    // Amplitude(s) for diagram number 785
    FFV1_0( w_fp[196], w_fp[229], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[229], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[229], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[15] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 786 OF 15495 ***
    // Wavefunction(s) for diagram number 786
    // (none)
    // Amplitude(s) for diagram number 786
    FFV1_0( w_fp[200], w_fp[2], w_fp[151], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[200], w_fp[2], w_fp[152], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[200], w_fp[2], w_fp[153], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 787 OF 15495 ***
    // Wavefunction(s) for diagram number 787
    // (none)
    // Amplitude(s) for diagram number 787
    FFV1_0( w_fp[196], w_fp[2], w_fp[67], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[154], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[11] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[17] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[155], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[15] -= amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 788 OF 15495 ***
    // Wavefunction(s) for diagram number 788
    // (none)
    // Amplitude(s) for diagram number 788
    FFV1_0( w_fp[218], w_fp[236], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 789 OF 15495 ***
    // Wavefunction(s) for diagram number 789
    // (none)
    // Amplitude(s) for diagram number 789
    FFV1_0( w_fp[218], w_fp[237], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 790 OF 15495 ***
    // Wavefunction(s) for diagram number 790
    FFV1_1( w_fp[229], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[155] );
    // Amplitude(s) for diagram number 790
    FFV1_0( w_fp[170], w_fp[155], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 791 OF 15495 ***
    // Wavefunction(s) for diagram number 791
    // (none)
    // Amplitude(s) for diagram number 791
    FFV1_0( w_fp[170], w_fp[237], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 792 OF 15495 ***
    // Wavefunction(s) for diagram number 792
    // (none)
    // Amplitude(s) for diagram number 792
    FFV1_0( w_fp[171], w_fp[155], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 793 OF 15495 ***
    // Wavefunction(s) for diagram number 793
    // (none)
    // Amplitude(s) for diagram number 793
    FFV1_0( w_fp[171], w_fp[236], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 794 OF 15495 ***
    // Wavefunction(s) for diagram number 794
    FFV1P0_3( w_fp[168], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[154] );
    // Amplitude(s) for diagram number 794
    VVVV1_0( w_fp[10], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
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
    VVVV3_0( w_fp[10], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0( w_fp[10], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 795 OF 15495 ***
    // Wavefunction(s) for diagram number 795
    // (none)
    // Amplitude(s) for diagram number 795
    VVV1_0( w_fp[154], w_fp[7], w_fp[12], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 796 OF 15495 ***
    // Wavefunction(s) for diagram number 796
    // (none)
    // Amplitude(s) for diagram number 796
    VVV1_0( w_fp[154], w_fp[6], w_fp[13], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 797 OF 15495 ***
    // Wavefunction(s) for diagram number 797
    FFV1_1( w_fp[2], w_fp[10], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[67] );
    // Amplitude(s) for diagram number 797
    FFV1_0( w_fp[170], w_fp[67], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];

    // *** DIAGRAM 798 OF 15495 ***
    // Wavefunction(s) for diagram number 798
    // (none)
    // Amplitude(s) for diagram number 798
    FFV1_0( w_fp[170], w_fp[2], w_fp[13], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 799 OF 15495 ***
    // Wavefunction(s) for diagram number 799
    // (none)
    // Amplitude(s) for diagram number 799
    FFV1_0( w_fp[171], w_fp[67], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];

    // *** DIAGRAM 800 OF 15495 ***
    // Wavefunction(s) for diagram number 800
    // (none)
    // Amplitude(s) for diagram number 800
    FFV1_0( w_fp[171], w_fp[2], w_fp[12], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 67 );
    storeWf( wfs, w_cx, nevt, 122 );
    storeWf( wfs, w_cx, nevt, 123 );
    storeWf( wfs, w_cx, nevt, 128 );
    storeWf( wfs, w_cx, nevt, 129 );
    storeWf( wfs, w_cx, nevt, 154 );
    storeWf( wfs, w_cx, nevt, 155 );
    storeWf( wfs, w_cx, nevt, 229 );
    storeWf( wfs, w_cx, nevt, 236 );
    storeWf( wfs, w_cx, nevt, 237 );
    storeWf( wfs, w_cx, nevt, 238 );
    storeWf( wfs, w_cx, nevt, 239 );
    storeWf( wfs, w_cx, nevt, 240 );
    storeWf( wfs, w_cx, nevt, 241 );
    storeWf( wfs, w_cx, nevt, 242 );
    storeWf( wfs, w_cx, nevt, 243 );
    storeWf( wfs, w_cx, nevt, 244 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
