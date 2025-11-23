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
  diagramgroup481( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 515 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 581 );
#endif
#endif

    // *** DIAGRAM 4801 OF 15495 ***
    // Wavefunction(s) for diagram number 4801
    // (none)
    // Amplitude(s) for diagram number 4801
    VVV1_0( w_fp[1], w_fp[144], w_fp[497], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 4802 OF 15495 ***
    // Wavefunction(s) for diagram number 4802
    // (none)
    // Amplitude(s) for diagram number 4802
    FFV1_0( w_fp[180], w_fp[2], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4803 OF 15495 ***
    // Wavefunction(s) for diagram number 4803
    // (none)
    // Amplitude(s) for diagram number 4803
    FFV1_0( w_fp[180], w_fp[536], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 4804 OF 15495 ***
    // Wavefunction(s) for diagram number 4804
    // (none)
    // Amplitude(s) for diagram number 4804
    VVV1_0( w_fp[576], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    VVV1_0( w_fp[577], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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
    VVV1_0( w_fp[578], w_fp[144], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 4805 OF 15495 ***
    // Wavefunction(s) for diagram number 4805
    // (none)
    // Amplitude(s) for diagram number 4805
    FFV1_0( w_fp[181], w_fp[2], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[181], w_fp[2], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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
    FFV1_0( w_fp[181], w_fp[2], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4806 OF 15495 ***
    // Wavefunction(s) for diagram number 4806
    // (none)
    // Amplitude(s) for diagram number 4806
    VVV1_0( w_fp[579], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    VVV1_0( w_fp[580], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    VVV1_0( w_fp[581], w_fp[144], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];

    // *** DIAGRAM 4807 OF 15495 ***
    // Wavefunction(s) for diagram number 4807
    // (none)
    // Amplitude(s) for diagram number 4807
    FFV1_0( w_fp[180], w_fp[2], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[180], w_fp[2], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[180], w_fp[2], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4808 OF 15495 ***
    // Wavefunction(s) for diagram number 4808
    // (none)
    // Amplitude(s) for diagram number 4808
    FFV1_0( w_fp[252], w_fp[2], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[252], w_fp[2], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[252], w_fp[2], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 4809 OF 15495 ***
    // Wavefunction(s) for diagram number 4809
    // (none)
    // Amplitude(s) for diagram number 4809
    VVV1_0( w_fp[515], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    VVV1_0( w_fp[454], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    VVV1_0( w_fp[438], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
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

    // *** DIAGRAM 4810 OF 15495 ***
    // Wavefunction(s) for diagram number 4810
    // (none)
    // Amplitude(s) for diagram number 4810
    VVV1_0( w_fp[546], w_fp[144], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------

#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup482( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 96 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 546 );
#endif
#endif

    // *** DIAGRAM 4811 OF 15495 ***
    // Wavefunction(s) for diagram number 4811
    // (none)
    // Amplitude(s) for diagram number 4811
    FFV1_0( w_fp[179], w_fp[244], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4812 OF 15495 ***
    // Wavefunction(s) for diagram number 4812
    // (none)
    // Amplitude(s) for diagram number 4812
    FFV1_0( w_fp[96], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4813 OF 15495 ***
    // Wavefunction(s) for diagram number 4813
    // (none)
    // Amplitude(s) for diagram number 4813
    FFV1_0( w_fp[252], w_fp[530], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];

    // *** DIAGRAM 4814 OF 15495 ***
    // Wavefunction(s) for diagram number 4814
    // (none)
    // Amplitude(s) for diagram number 4814
    FFV1_0( w_fp[179], w_fp[530], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4815 OF 15495 ***
    // Wavefunction(s) for diagram number 4815
    // (none)
    // Amplitude(s) for diagram number 4815
    FFV1_0( w_fp[96], w_fp[530], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];

    // *** DIAGRAM 4816 OF 15495 ***
    // Wavefunction(s) for diagram number 4816
    // (none)
    // Amplitude(s) for diagram number 4816
    FFV1_0( w_fp[516], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4817 OF 15495 ***
    // Wavefunction(s) for diagram number 4817
    // (none)
    // Amplitude(s) for diagram number 4817
    FFV1_0( w_fp[516], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];

    // *** DIAGRAM 4818 OF 15495 ***
    // Wavefunction(s) for diagram number 4818
    // (none)
    // Amplitude(s) for diagram number 4818
    FFV1_0( w_fp[252], w_fp[2], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4819 OF 15495 ***
    // Wavefunction(s) for diagram number 4819
    // (none)
    // Amplitude(s) for diagram number 4819
    VVV1_0( w_fp[518], w_fp[1], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 4820 OF 15495 ***
    // Wavefunction(s) for diagram number 4820
    // (none)
    // Amplitude(s) for diagram number 4820
    FFV1_0( w_fp[252], w_fp[244], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup483( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 589 );
    retrieveWf( wfs, w_cx, nevt, 590 );
#endif
#endif

    // *** DIAGRAM 4821 OF 15495 ***
    // Wavefunction(s) for diagram number 4821
    // (none)
    // Amplitude(s) for diagram number 4821
    VVV1_0( w_fp[435], w_fp[359], w_fp[144], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 4822 OF 15495 ***
    // Wavefunction(s) for diagram number 4822
    // (none)
    // Amplitude(s) for diagram number 4822
    FFV1_0( w_fp[179], w_fp[2], w_fp[588], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[204] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    FFV1_0( w_fp[179], w_fp[2], w_fp[590], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 4823 OF 15495 ***
    // Wavefunction(s) for diagram number 4823
    // (none)
    // Amplitude(s) for diagram number 4823
    FFV1_0( w_fp[529], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[574] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4824 OF 15495 ***
    // Wavefunction(s) for diagram number 4824
    // (none)
    // Amplitude(s) for diagram number 4824
    FFV1_0( w_fp[3], w_fp[244], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];

    // *** DIAGRAM 4825 OF 15495 ***
    // Wavefunction(s) for diagram number 4825
    // (none)
    // Amplitude(s) for diagram number 4825
    FFV1_0( w_fp[184], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4826 OF 15495 ***
    // Wavefunction(s) for diagram number 4826
    // (none)
    // Amplitude(s) for diagram number 4826
    FFV1_0( w_fp[184], w_fp[2], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4827 OF 15495 ***
    // Wavefunction(s) for diagram number 4827
    // (none)
    // Amplitude(s) for diagram number 4827
    FFV1_0( w_fp[3], w_fp[453], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 4828 OF 15495 ***
    // Wavefunction(s) for diagram number 4828
    // (none)
    // Amplitude(s) for diagram number 4828
    FFV1_0( w_fp[529], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 4829 OF 15495 ***
    // Wavefunction(s) for diagram number 4829
    // (none)
    // Amplitude(s) for diagram number 4829
    VVV1_0( w_fp[359], w_fp[7], w_fp[538], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];

    // *** DIAGRAM 4830 OF 15495 ***
    // Wavefunction(s) for diagram number 4830
    // (none)
    // Amplitude(s) for diagram number 4830
    FFV1_0( w_fp[3], w_fp[532], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[45] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[261] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[263] += cxtype( 0, 1 ) * amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup484( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 374 );
    retrieveWf( wfs, w_cx, nevt, 375 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 545 );
#endif
#endif

    // *** DIAGRAM 4831 OF 15495 ***
    // Wavefunction(s) for diagram number 4831
    // (none)
    // Amplitude(s) for diagram number 4831
    FFV1_0( w_fp[184], w_fp[527], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 4832 OF 15495 ***
    // Wavefunction(s) for diagram number 4832
    // (none)
    // Amplitude(s) for diagram number 4832
    FFV1_0( w_fp[184], w_fp[532], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];

    // *** DIAGRAM 4833 OF 15495 ***
    // Wavefunction(s) for diagram number 4833
    // (none)
    // Amplitude(s) for diagram number 4833
    FFV1_0( w_fp[3], w_fp[527], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4834 OF 15495 ***
    // Wavefunction(s) for diagram number 4834
    // (none)
    // Amplitude(s) for diagram number 4834
    VVV1_0( w_fp[1], w_fp[116], w_fp[538], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];

    // *** DIAGRAM 4835 OF 15495 ***
    // Wavefunction(s) for diagram number 4835
    // (none)
    // Amplitude(s) for diagram number 4835
    FFV1_0( w_fp[3], w_fp[530], w_fp[79], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[530], w_fp[374], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[530], w_fp[375], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];

    // *** DIAGRAM 4836 OF 15495 ***
    // Wavefunction(s) for diagram number 4836
    // (none)
    // Amplitude(s) for diagram number 4836
    VVV1_0( w_fp[359], w_fp[7], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 4837 OF 15495 ***
    // Wavefunction(s) for diagram number 4837
    // (none)
    // Amplitude(s) for diagram number 4837
    FFV1_0( w_fp[442], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[209] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4838 OF 15495 ***
    // Wavefunction(s) for diagram number 4838
    // (none)
    // Amplitude(s) for diagram number 4838
    FFV1_0( w_fp[528], w_fp[244], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];

    // *** DIAGRAM 4839 OF 15495 ***
    // Wavefunction(s) for diagram number 4839
    // (none)
    // Amplitude(s) for diagram number 4839
    FFV1_0( w_fp[442], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[442] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];

    // *** DIAGRAM 4840 OF 15495 ***
    // Wavefunction(s) for diagram number 4840
    // (none)
    // Amplitude(s) for diagram number 4840
    FFV1_0( w_fp[528], w_fp[2], w_fp[116], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup485( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 79 );
    retrieveWf( wfs, w_cx, nevt, 116 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 374 );
    retrieveWf( wfs, w_cx, nevt, 375 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
#endif
#endif

    // *** DIAGRAM 4841 OF 15495 ***
    // Wavefunction(s) for diagram number 4841
    // (none)
    // Amplitude(s) for diagram number 4841
    VVV1_0( w_fp[1], w_fp[116], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 4842 OF 15495 ***
    // Wavefunction(s) for diagram number 4842
    // (none)
    // Amplitude(s) for diagram number 4842
    FFV1_0( w_fp[479], w_fp[2], w_fp[79], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    FFV1_0( w_fp[479], w_fp[2], w_fp[374], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    FFV1_0( w_fp[479], w_fp[2], w_fp[375], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 4843 OF 15495 ***
    // Wavefunction(s) for diagram number 4843
    // (none)
    // Amplitude(s) for diagram number 4843
    FFV1_0( w_fp[3], w_fp[537], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[45] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[261] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[609] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[651] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];

    // *** DIAGRAM 4844 OF 15495 ***
    // Wavefunction(s) for diagram number 4844
    // (none)
    // Amplitude(s) for diagram number 4844
    FFV1_0( w_fp[486], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[209] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];

    // *** DIAGRAM 4845 OF 15495 ***
    // Wavefunction(s) for diagram number 4845
    // (none)
    // Amplitude(s) for diagram number 4845
    FFV1_0( w_fp[3], w_fp[244], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[555] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[574] += amp_sv[0];

    // *** DIAGRAM 4846 OF 15495 ***
    // Wavefunction(s) for diagram number 4846
    // (none)
    // Amplitude(s) for diagram number 4846
    FFV1_0( w_fp[486], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[563] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4847 OF 15495 ***
    // Wavefunction(s) for diagram number 4847
    // (none)
    // Amplitude(s) for diagram number 4847
    FFV1_0( w_fp[184], w_fp[2], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[43] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4848 OF 15495 ***
    // Wavefunction(s) for diagram number 4848
    // (none)
    // Amplitude(s) for diagram number 4848
    FFV1_0( w_fp[184], w_fp[537], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4849 OF 15495 ***
    // Wavefunction(s) for diagram number 4849
    // (none)
    // Amplitude(s) for diagram number 4849
    FFV1_0( w_fp[3], w_fp[244], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[244], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[555] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[563] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[244], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];

    // *** DIAGRAM 4850 OF 15495 ***
    // Wavefunction(s) for diagram number 4850
    // (none)
    // Amplitude(s) for diagram number 4850
    FFV1_0( w_fp[184], w_fp[2], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    FFV1_0( w_fp[184], w_fp[2], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    FFV1_0( w_fp[184], w_fp[2], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup486( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 561 );
#endif
#endif

    // *** DIAGRAM 4851 OF 15495 ***
    // Wavefunction(s) for diagram number 4851
    // (none)
    // Amplitude(s) for diagram number 4851
    FFV1_0( w_fp[529], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[475] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[694] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4852 OF 15495 ***
    // Wavefunction(s) for diagram number 4852
    // (none)
    // Amplitude(s) for diagram number 4852
    FFV1_0( w_fp[3], w_fp[122], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 4853 OF 15495 ***
    // Wavefunction(s) for diagram number 4853
    // (none)
    // Amplitude(s) for diagram number 4853
    FFV1_0( w_fp[186], w_fp[453], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4854 OF 15495 ***
    // Wavefunction(s) for diagram number 4854
    // (none)
    // Amplitude(s) for diagram number 4854
    FFV1_0( w_fp[186], w_fp[2], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 4855 OF 15495 ***
    // Wavefunction(s) for diagram number 4855
    // (none)
    // Amplitude(s) for diagram number 4855
    FFV1_0( w_fp[3], w_fp[453], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];

    // *** DIAGRAM 4856 OF 15495 ***
    // Wavefunction(s) for diagram number 4856
    // (none)
    // Amplitude(s) for diagram number 4856
    FFV1_0( w_fp[529], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[574] -= amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[598] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];

    // *** DIAGRAM 4857 OF 15495 ***
    // Wavefunction(s) for diagram number 4857
    // (none)
    // Amplitude(s) for diagram number 4857
    VVV1_0( w_fp[360], w_fp[6], w_fp[538], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];

    // *** DIAGRAM 4858 OF 15495 ***
    // Wavefunction(s) for diagram number 4858
    // (none)
    // Amplitude(s) for diagram number 4858
    FFV1_0( w_fp[3], w_fp[531], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[255] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[257] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4859 OF 15495 ***
    // Wavefunction(s) for diagram number 4859
    // (none)
    // Amplitude(s) for diagram number 4859
    FFV1_0( w_fp[186], w_fp[527], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];

    // *** DIAGRAM 4860 OF 15495 ***
    // Wavefunction(s) for diagram number 4860
    // (none)
    // Amplitude(s) for diagram number 4860
    FFV1_0( w_fp[186], w_fp[531], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup487( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 46 );
    retrieveWf( wfs, w_cx, nevt, 81 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 125 );
    retrieveWf( wfs, w_cx, nevt, 358 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 545 );
#endif
#endif

    // *** DIAGRAM 4861 OF 15495 ***
    // Wavefunction(s) for diagram number 4861
    // (none)
    // Amplitude(s) for diagram number 4861
    FFV1_0( w_fp[3], w_fp[527], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4862 OF 15495 ***
    // Wavefunction(s) for diagram number 4862
    // (none)
    // Amplitude(s) for diagram number 4862
    VVV1_0( w_fp[1], w_fp[125], w_fp[538], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];

    // *** DIAGRAM 4863 OF 15495 ***
    // Wavefunction(s) for diagram number 4863
    // (none)
    // Amplitude(s) for diagram number 4863
    FFV1_0( w_fp[3], w_fp[530], w_fp[358], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[255] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[530], w_fp[81], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[45] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[261] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[530], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];

    // *** DIAGRAM 4864 OF 15495 ***
    // Wavefunction(s) for diagram number 4864
    // (none)
    // Amplitude(s) for diagram number 4864
    VVV1_0( w_fp[360], w_fp[6], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[521] -= amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 4865 OF 15495 ***
    // Wavefunction(s) for diagram number 4865
    // (none)
    // Amplitude(s) for diagram number 4865
    FFV1_0( w_fp[505], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[232] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4866 OF 15495 ***
    // Wavefunction(s) for diagram number 4866
    // (none)
    // Amplitude(s) for diagram number 4866
    FFV1_0( w_fp[528], w_fp[122], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];

    // *** DIAGRAM 4867 OF 15495 ***
    // Wavefunction(s) for diagram number 4867
    // (none)
    // Amplitude(s) for diagram number 4867
    FFV1_0( w_fp[505], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[466] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];

    // *** DIAGRAM 4868 OF 15495 ***
    // Wavefunction(s) for diagram number 4868
    // (none)
    // Amplitude(s) for diagram number 4868
    FFV1_0( w_fp[528], w_fp[2], w_fp[125], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4869 OF 15495 ***
    // Wavefunction(s) for diagram number 4869
    // (none)
    // Amplitude(s) for diagram number 4869
    VVV1_0( w_fp[1], w_fp[125], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];

    // *** DIAGRAM 4870 OF 15495 ***
    // Wavefunction(s) for diagram number 4870
    // (none)
    // Amplitude(s) for diagram number 4870
    FFV1_0( w_fp[479], w_fp[2], w_fp[358], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    FFV1_0( w_fp[479], w_fp[2], w_fp[81], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[208] += amp_sv[0];
    jamp_sv[209] -= amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    FFV1_0( w_fp[479], w_fp[2], w_fp[46], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[527] -= amp_sv[0];
    jamp_sv[572] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup488( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 581 );
#endif
#endif

    // *** DIAGRAM 4871 OF 15495 ***
    // Wavefunction(s) for diagram number 4871
    // (none)
    // Amplitude(s) for diagram number 4871
    FFV1_0( w_fp[3], w_fp[536], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[255] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[531] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];

    // *** DIAGRAM 4872 OF 15495 ***
    // Wavefunction(s) for diagram number 4872
    // (none)
    // Amplitude(s) for diagram number 4872
    FFV1_0( w_fp[450], w_fp[2], w_fp[360], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];

    // *** DIAGRAM 4873 OF 15495 ***
    // Wavefunction(s) for diagram number 4873
    // (none)
    // Amplitude(s) for diagram number 4873
    FFV1_0( w_fp[3], w_fp[122], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[694] += amp_sv[0];

    // *** DIAGRAM 4874 OF 15495 ***
    // Wavefunction(s) for diagram number 4874
    // (none)
    // Amplitude(s) for diagram number 4874
    FFV1_0( w_fp[450], w_fp[122], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[683] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4875 OF 15495 ***
    // Wavefunction(s) for diagram number 4875
    // (none)
    // Amplitude(s) for diagram number 4875
    FFV1_0( w_fp[186], w_fp[2], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 4876 OF 15495 ***
    // Wavefunction(s) for diagram number 4876
    // (none)
    // Amplitude(s) for diagram number 4876
    FFV1_0( w_fp[186], w_fp[536], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4877 OF 15495 ***
    // Wavefunction(s) for diagram number 4877
    // (none)
    // Amplitude(s) for diagram number 4877
    FFV1_0( w_fp[3], w_fp[122], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[692] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[122], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[683] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[122], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[692] -= amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[694] -= amp_sv[0];

    // *** DIAGRAM 4878 OF 15495 ***
    // Wavefunction(s) for diagram number 4878
    // (none)
    // Amplitude(s) for diagram number 4878
    FFV1_0( w_fp[186], w_fp[2], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[242] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    FFV1_0( w_fp[186], w_fp[2], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    FFV1_0( w_fp[186], w_fp[2], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];

    // *** DIAGRAM 4879 OF 15495 ***
    // Wavefunction(s) for diagram number 4879
    // (none)
    // Amplitude(s) for diagram number 4879
    FFV1_0( w_fp[529], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[598] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[718] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4880 OF 15495 ***
    // Wavefunction(s) for diagram number 4880
    // (none)
    // Amplitude(s) for diagram number 4880
    FFV1_0( w_fp[3], w_fp[128], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup489( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 549 );
#endif
#endif

    // *** DIAGRAM 4881 OF 15495 ***
    // Wavefunction(s) for diagram number 4881
    // (none)
    // Amplitude(s) for diagram number 4881
    FFV1_0( w_fp[188], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4882 OF 15495 ***
    // Wavefunction(s) for diagram number 4882
    // (none)
    // Amplitude(s) for diagram number 4882
    FFV1_0( w_fp[188], w_fp[2], w_fp[549], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];

    // *** DIAGRAM 4883 OF 15495 ***
    // Wavefunction(s) for diagram number 4883
    // (none)
    // Amplitude(s) for diagram number 4883
    FFV1_0( w_fp[3], w_fp[453], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 4884 OF 15495 ***
    // Wavefunction(s) for diagram number 4884
    // (none)
    // Amplitude(s) for diagram number 4884
    FFV1_0( w_fp[529], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[598] -= amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[718] += amp_sv[0];

    // *** DIAGRAM 4885 OF 15495 ***
    // Wavefunction(s) for diagram number 4885
    // (none)
    // Amplitude(s) for diagram number 4885
    VVV1_0( w_fp[361], w_fp[5], w_fp[538], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[35] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[249] += amp_sv[0];
    jamp_sv[251] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];

    // *** DIAGRAM 4886 OF 15495 ***
    // Wavefunction(s) for diagram number 4886
    // (none)
    // Amplitude(s) for diagram number 4886
    FFV1_0( w_fp[3], w_fp[533], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[35] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[249] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[251] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4887 OF 15495 ***
    // Wavefunction(s) for diagram number 4887
    // (none)
    // Amplitude(s) for diagram number 4887
    FFV1_0( w_fp[188], w_fp[527], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];

    // *** DIAGRAM 4888 OF 15495 ***
    // Wavefunction(s) for diagram number 4888
    // (none)
    // Amplitude(s) for diagram number 4888
    FFV1_0( w_fp[188], w_fp[533], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];

    // *** DIAGRAM 4889 OF 15495 ***
    // Wavefunction(s) for diagram number 4889
    // (none)
    // Amplitude(s) for diagram number 4889
    FFV1_0( w_fp[3], w_fp[527], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4890 OF 15495 ***
    // Wavefunction(s) for diagram number 4890
    // (none)
    // Amplitude(s) for diagram number 4890
    VVV1_0( w_fp[1], w_fp[131], w_fp[538], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup490( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 131 );
    retrieveWf( wfs, w_cx, nevt, 249 );
    retrieveWf( wfs, w_cx, nevt, 257 );
    retrieveWf( wfs, w_cx, nevt, 354 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 545 );
#endif
#endif

    // *** DIAGRAM 4891 OF 15495 ***
    // Wavefunction(s) for diagram number 4891
    // (none)
    // Amplitude(s) for diagram number 4891
    FFV1_0( w_fp[3], w_fp[530], w_fp[257], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[530], w_fp[249], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[530], w_fp[354], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];

    // *** DIAGRAM 4892 OF 15495 ***
    // Wavefunction(s) for diagram number 4892
    // (none)
    // Amplitude(s) for diagram number 4892
    VVV1_0( w_fp[361], w_fp[5], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 4893 OF 15495 ***
    // Wavefunction(s) for diagram number 4893
    // (none)
    // Amplitude(s) for diagram number 4893
    FFV1_0( w_fp[447], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[587] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[707] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4894 OF 15495 ***
    // Wavefunction(s) for diagram number 4894
    // (none)
    // Amplitude(s) for diagram number 4894
    FFV1_0( w_fp[528], w_fp[128], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 4895 OF 15495 ***
    // Wavefunction(s) for diagram number 4895
    // (none)
    // Amplitude(s) for diagram number 4895
    FFV1_0( w_fp[447], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[586] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[707] += amp_sv[0];

    // *** DIAGRAM 4896 OF 15495 ***
    // Wavefunction(s) for diagram number 4896
    // (none)
    // Amplitude(s) for diagram number 4896
    FFV1_0( w_fp[528], w_fp[2], w_fp[131], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[453] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[477] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4897 OF 15495 ***
    // Wavefunction(s) for diagram number 4897
    // (none)
    // Amplitude(s) for diagram number 4897
    VVV1_0( w_fp[1], w_fp[131], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];

    // *** DIAGRAM 4898 OF 15495 ***
    // Wavefunction(s) for diagram number 4898
    // (none)
    // Amplitude(s) for diagram number 4898
    FFV1_0( w_fp[479], w_fp[2], w_fp[257], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[716] += amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    FFV1_0( w_fp[479], w_fp[2], w_fp[249], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[453] -= amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    jamp_sv[477] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    FFV1_0( w_fp[479], w_fp[2], w_fp[354], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 4899 OF 15495 ***
    // Wavefunction(s) for diagram number 4899
    // (none)
    // Amplitude(s) for diagram number 4899
    FFV1_0( w_fp[3], w_fp[535], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[35] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[249] -= amp_sv[0];
    jamp_sv[251] += amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];

    // *** DIAGRAM 4900 OF 15495 ***
    // Wavefunction(s) for diagram number 4900
    // (none)
    // Amplitude(s) for diagram number 4900
    FFV1_0( w_fp[495], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];

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
    // (none)
#endif
#endif
  }

  //--------------------------------------------------------------------------

}
