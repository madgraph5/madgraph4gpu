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
  diagramgroup711( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 580 );
#endif
#endif

    // *** DIAGRAM 7101 OF 15495 ***
    // Wavefunction(s) for diagram number 7101
    // (none)
    // Amplitude(s) for diagram number 7101
    FFV1_0( w_fp[521], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7102 OF 15495 ***
    // Wavefunction(s) for diagram number 7102
    // (none)
    // Amplitude(s) for diagram number 7102
    FFV1_0( w_fp[562], w_fp[225], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7103 OF 15495 ***
    // Wavefunction(s) for diagram number 7103
    // (none)
    // Amplitude(s) for diagram number 7103
    FFV1_0( w_fp[521], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7104 OF 15495 ***
    // Wavefunction(s) for diagram number 7104
    // (none)
    // Amplitude(s) for diagram number 7104
    FFV1_0( w_fp[562], w_fp[226], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7105 OF 15495 ***
    // Wavefunction(s) for diagram number 7105
    // (none)
    // Amplitude(s) for diagram number 7105
    FFV1_0( w_fp[524], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7106 OF 15495 ***
    // Wavefunction(s) for diagram number 7106
    // (none)
    // Amplitude(s) for diagram number 7106
    FFV1_0( w_fp[577], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[637] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];

    // *** DIAGRAM 7107 OF 15495 ***
    // Wavefunction(s) for diagram number 7107
    // (none)
    // Amplitude(s) for diagram number 7107
    FFV1_0( w_fp[3], w_fp[473], w_fp[580], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7108 OF 15495 ***
    // Wavefunction(s) for diagram number 7108
    // (none)
    // Amplitude(s) for diagram number 7108
    VVVV1_0( w_fp[523], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    VVVV3_0( w_fp[523], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    VVVV4_0( w_fp[523], w_fp[1], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 7109 OF 15495 ***
    // Wavefunction(s) for diagram number 7109
    // (none)
    // Amplitude(s) for diagram number 7109
    VVV1_0( w_fp[228], w_fp[5], w_fp[578], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];

    // *** DIAGRAM 7110 OF 15495 ***
    // Wavefunction(s) for diagram number 7110
    // (none)
    // Amplitude(s) for diagram number 7110
    VVV1_0( w_fp[1], w_fp[228], w_fp[580], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

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
  diagramgroup712( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 582 );
#endif
#endif

    // *** DIAGRAM 7111 OF 15495 ***
    // Wavefunction(s) for diagram number 7111
    // (none)
    // Amplitude(s) for diagram number 7111
    FFV1_0( w_fp[3], w_fp[226], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7112 OF 15495 ***
    // Wavefunction(s) for diagram number 7112
    // (none)
    // Amplitude(s) for diagram number 7112
    FFV1_0( w_fp[577], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[679] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];

    // *** DIAGRAM 7113 OF 15495 ***
    // Wavefunction(s) for diagram number 7113
    // (none)
    // Amplitude(s) for diagram number 7113
    FFV1_0( w_fp[514], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[631] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];

    // *** DIAGRAM 7114 OF 15495 ***
    // Wavefunction(s) for diagram number 7114
    // (none)
    // Amplitude(s) for diagram number 7114
    FFV1_0( w_fp[3], w_fp[473], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7115 OF 15495 ***
    // Wavefunction(s) for diagram number 7115
    // (none)
    // Amplitude(s) for diagram number 7115
    VVVV1_0( w_fp[474], w_fp[1], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[1], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[1], w_fp[228], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 7116 OF 15495 ***
    // Wavefunction(s) for diagram number 7116
    // (none)
    // Amplitude(s) for diagram number 7116
    VVV1_0( w_fp[228], w_fp[4], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 7117 OF 15495 ***
    // Wavefunction(s) for diagram number 7117
    // (none)
    // Amplitude(s) for diagram number 7117
    VVV1_0( w_fp[1], w_fp[228], w_fp[479], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 7118 OF 15495 ***
    // Wavefunction(s) for diagram number 7118
    // (none)
    // Amplitude(s) for diagram number 7118
    FFV1_0( w_fp[3], w_fp[225], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[653] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7119 OF 15495 ***
    // Wavefunction(s) for diagram number 7119
    // (none)
    // Amplitude(s) for diagram number 7119
    FFV1_0( w_fp[514], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[655] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 7120 OF 15495 ***
    // Wavefunction(s) for diagram number 7120
    // (none)
    // Amplitude(s) for diagram number 7120
    VVV1_0( w_fp[516], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    VVV1_0( w_fp[537], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[681] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    VVV1_0( w_fp[582], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];

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
  diagramgroup713( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 230 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
#endif
#endif

    // *** DIAGRAM 7121 OF 15495 ***
    // Wavefunction(s) for diagram number 7121
    // (none)
    // Amplitude(s) for diagram number 7121
    FFV1_0( w_fp[3], w_fp[226], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7122 OF 15495 ***
    // Wavefunction(s) for diagram number 7122
    // (none)
    // Amplitude(s) for diagram number 7122
    VVV1_0( w_fp[525], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    VVV1_0( w_fp[526], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    VVV1_0( w_fp[45], w_fp[228], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7123 OF 15495 ***
    // Wavefunction(s) for diagram number 7123
    // (none)
    // Amplitude(s) for diagram number 7123
    FFV1_0( w_fp[3], w_fp[225], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[225], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[225], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7124 OF 15495 ***
    // Wavefunction(s) for diagram number 7124
    // (none)
    // Amplitude(s) for diagram number 7124
    FFV1_0( w_fp[3], w_fp[473], w_fp[588], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[586], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[473], w_fp[451], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7125 OF 15495 ***
    // Wavefunction(s) for diagram number 7125
    // (none)
    // Amplitude(s) for diagram number 7125
    VVV1_0( w_fp[588], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0( w_fp[586], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[623] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    VVV1_0( w_fp[451], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];

    // *** DIAGRAM 7126 OF 15495 ***
    // Wavefunction(s) for diagram number 7126
    // (none)
    // Amplitude(s) for diagram number 7126
    VVV1_0( w_fp[450], w_fp[230], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7127 OF 15495 ***
    // Wavefunction(s) for diagram number 7127
    // (none)
    // Amplitude(s) for diagram number 7127
    FFV1_0( w_fp[196], w_fp[226], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];

    // *** DIAGRAM 7128 OF 15495 ***
    // Wavefunction(s) for diagram number 7128
    // (none)
    // Amplitude(s) for diagram number 7128
    FFV1_0( w_fp[216], w_fp[215], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];

    // *** DIAGRAM 7129 OF 15495 ***
    // Wavefunction(s) for diagram number 7129
    // (none)
    // Amplitude(s) for diagram number 7129
    FFV1_0( w_fp[355], w_fp[435], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7130 OF 15495 ***
    // Wavefunction(s) for diagram number 7130
    // (none)
    // Amplitude(s) for diagram number 7130
    FFV1_0( w_fp[216], w_fp[435], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup714( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 230 );
    retrieveWf( wfs, w_cx, nevt, 231 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 530 );
#endif
#endif

    // *** DIAGRAM 7131 OF 15495 ***
    // Wavefunction(s) for diagram number 7131
    // (none)
    // Amplitude(s) for diagram number 7131
    FFV1_0( w_fp[529], w_fp[473], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7132 OF 15495 ***
    // Wavefunction(s) for diagram number 7132
    // (none)
    // Amplitude(s) for diagram number 7132
    FFV1_0( w_fp[529], w_fp[226], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7133 OF 15495 ***
    // Wavefunction(s) for diagram number 7133
    // (none)
    // Amplitude(s) for diagram number 7133
    FFV1_0( w_fp[196], w_fp[473], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[629] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];

    // *** DIAGRAM 7134 OF 15495 ***
    // Wavefunction(s) for diagram number 7134
    // (none)
    // Amplitude(s) for diagram number 7134
    FFV1_0( w_fp[355], w_fp[215], w_fp[474], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 7135 OF 15495 ***
    // Wavefunction(s) for diagram number 7135
    // (none)
    // Amplitude(s) for diagram number 7135
    VVV1_0( w_fp[474], w_fp[1], w_fp[230], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7136 OF 15495 ***
    // Wavefunction(s) for diagram number 7136
    // (none)
    // Amplitude(s) for diagram number 7136
    FFV1_0( w_fp[216], w_fp[473], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7137 OF 15495 ***
    // Wavefunction(s) for diagram number 7137
    // (none)
    // Amplitude(s) for diagram number 7137
    FFV1_0( w_fp[355], w_fp[226], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7138 OF 15495 ***
    // Wavefunction(s) for diagram number 7138
    // (none)
    // Amplitude(s) for diagram number 7138
    FFV1_0( w_fp[196], w_fp[215], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[215], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[640] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[215], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7139 OF 15495 ***
    // Wavefunction(s) for diagram number 7139
    // (none)
    // Amplitude(s) for diagram number 7139
    VVV1_0( w_fp[450], w_fp[231], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7140 OF 15495 ***
    // Wavefunction(s) for diagram number 7140
    // (none)
    // Amplitude(s) for diagram number 7140
    FFV1_0( w_fp[168], w_fp[225], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[652] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];

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
  diagramgroup715( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 231 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 555 );
#endif
#endif

    // *** DIAGRAM 7141 OF 15495 ***
    // Wavefunction(s) for diagram number 7141
    // (none)
    // Amplitude(s) for diagram number 7141
    FFV1_0( w_fp[218], w_fp[215], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];

    // *** DIAGRAM 7142 OF 15495 ***
    // Wavefunction(s) for diagram number 7142
    // (none)
    // Amplitude(s) for diagram number 7142
    FFV1_0( w_fp[256], w_fp[435], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7143 OF 15495 ***
    // Wavefunction(s) for diagram number 7143
    // (none)
    // Amplitude(s) for diagram number 7143
    FFV1_0( w_fp[218], w_fp[435], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7144 OF 15495 ***
    // Wavefunction(s) for diagram number 7144
    // (none)
    // Amplitude(s) for diagram number 7144
    FFV1_0( w_fp[555], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7145 OF 15495 ***
    // Wavefunction(s) for diagram number 7145
    // (none)
    // Amplitude(s) for diagram number 7145
    FFV1_0( w_fp[555], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7146 OF 15495 ***
    // Wavefunction(s) for diagram number 7146
    // (none)
    // Amplitude(s) for diagram number 7146
    FFV1_0( w_fp[168], w_fp[473], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];

    // *** DIAGRAM 7147 OF 15495 ***
    // Wavefunction(s) for diagram number 7147
    // (none)
    // Amplitude(s) for diagram number 7147
    FFV1_0( w_fp[256], w_fp[215], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];

    // *** DIAGRAM 7148 OF 15495 ***
    // Wavefunction(s) for diagram number 7148
    // (none)
    // Amplitude(s) for diagram number 7148
    VVV1_0( w_fp[523], w_fp[1], w_fp[231], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7149 OF 15495 ***
    // Wavefunction(s) for diagram number 7149
    // (none)
    // Amplitude(s) for diagram number 7149
    FFV1_0( w_fp[218], w_fp[473], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7150 OF 15495 ***
    // Wavefunction(s) for diagram number 7150
    // (none)
    // Amplitude(s) for diagram number 7150
    FFV1_0( w_fp[256], w_fp[225], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup716( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 233 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 582 );
#endif
#endif

    // *** DIAGRAM 7151 OF 15495 ***
    // Wavefunction(s) for diagram number 7151
    // (none)
    // Amplitude(s) for diagram number 7151
    FFV1_0( w_fp[168], w_fp[215], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[537], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7152 OF 15495 ***
    // Wavefunction(s) for diagram number 7152
    // (none)
    // Amplitude(s) for diagram number 7152
    VVV1_0( w_fp[450], w_fp[228], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7153 OF 15495 ***
    // Wavefunction(s) for diagram number 7153
    // (none)
    // Amplitude(s) for diagram number 7153
    FFV1_0( w_fp[3], w_fp[233], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[664] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[688] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7154 OF 15495 ***
    // Wavefunction(s) for diagram number 7154
    // (none)
    // Amplitude(s) for diagram number 7154
    FFV1_0( w_fp[221], w_fp[215], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7155 OF 15495 ***
    // Wavefunction(s) for diagram number 7155
    // (none)
    // Amplitude(s) for diagram number 7155
    FFV1_0( w_fp[3], w_fp[435], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[621] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7156 OF 15495 ***
    // Wavefunction(s) for diagram number 7156
    // (none)
    // Amplitude(s) for diagram number 7156
    FFV1_0( w_fp[221], w_fp[435], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7157 OF 15495 ***
    // Wavefunction(s) for diagram number 7157
    // (none)
    // Amplitude(s) for diagram number 7157
    FFV1_0( w_fp[532], w_fp[473], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];

    // *** DIAGRAM 7158 OF 15495 ***
    // Wavefunction(s) for diagram number 7158
    // (none)
    // Amplitude(s) for diagram number 7158
    FFV1_0( w_fp[532], w_fp[215], w_fp[254], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7159 OF 15495 ***
    // Wavefunction(s) for diagram number 7159
    // (none)
    // Amplitude(s) for diagram number 7159
    FFV1_0( w_fp[532], w_fp[233], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[662] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];

    // *** DIAGRAM 7160 OF 15495 ***
    // Wavefunction(s) for diagram number 7160
    // (none)
    // Amplitude(s) for diagram number 7160
    FFV1_0( w_fp[3], w_fp[473], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup717( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 579 );
#endif
#endif

    // *** DIAGRAM 7161 OF 15495 ***
    // Wavefunction(s) for diagram number 7161
    // (none)
    // Amplitude(s) for diagram number 7161
    VVV1_0( w_fp[518], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[661] += amp_sv[0];
    jamp_sv[664] -= amp_sv[0];
    jamp_sv[685] -= amp_sv[0];
    jamp_sv[688] += amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];

    // *** DIAGRAM 7162 OF 15495 ***
    // Wavefunction(s) for diagram number 7162
    // (none)
    // Amplitude(s) for diagram number 7162
    FFV1_0( w_fp[221], w_fp[473], w_fp[530], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];

    // *** DIAGRAM 7163 OF 15495 ***
    // Wavefunction(s) for diagram number 7163
    // (none)
    // Amplitude(s) for diagram number 7163
    VVV1_0( w_fp[530], w_fp[254], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 7164 OF 15495 ***
    // Wavefunction(s) for diagram number 7164
    // (none)
    // Amplitude(s) for diagram number 7164
    FFV1_0( w_fp[3], w_fp[215], w_fp[509], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[444], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[621] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[452], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[661] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[663] -= amp_sv[0];
    jamp_sv[664] += amp_sv[0];
    jamp_sv[685] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[687] += amp_sv[0];
    jamp_sv[688] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7165 OF 15495 ***
    // Wavefunction(s) for diagram number 7165
    // (none)
    // Amplitude(s) for diagram number 7165
    VVVV1_0( w_fp[450], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    VVVV3_0( w_fp[450], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    VVVV4_0( w_fp[450], w_fp[239], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7166 OF 15495 ***
    // Wavefunction(s) for diagram number 7166
    // (none)
    // Amplitude(s) for diagram number 7166
    VVV1_0( w_fp[239], w_fp[7], w_fp[447], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7167 OF 15495 ***
    // Wavefunction(s) for diagram number 7167
    // (none)
    // Amplitude(s) for diagram number 7167
    VVV1_0( w_fp[239], w_fp[5], w_fp[579], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7168 OF 15495 ***
    // Wavefunction(s) for diagram number 7168
    FFV1_1( w_fp[2], w_fp[450], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[442] );
    // Amplitude(s) for diagram number 7168
    FFV1_0( w_fp[216], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];

    // *** DIAGRAM 7169 OF 15495 ***
    // Wavefunction(s) for diagram number 7169
    // (none)
    // Amplitude(s) for diagram number 7169
    FFV1_0( w_fp[216], w_fp[2], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7170 OF 15495 ***
    // Wavefunction(s) for diagram number 7170
    // (none)
    // Amplitude(s) for diagram number 7170
    FFV1_0( w_fp[199], w_fp[442], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 442 );
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup718( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 447 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 590 );
    retrieveWf( wfs, w_cx, nevt, 594 );
#endif
#endif

    // *** DIAGRAM 7171 OF 15495 ***
    // Wavefunction(s) for diagram number 7171
    // (none)
    // Amplitude(s) for diagram number 7171
    FFV1_0( w_fp[199], w_fp[2], w_fp[447], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7172 OF 15495 ***
    // Wavefunction(s) for diagram number 7172
    // (none)
    // Amplitude(s) for diagram number 7172
    FFV1_0( w_fp[355], w_fp[590], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7173 OF 15495 ***
    // Wavefunction(s) for diagram number 7173
    // (none)
    // Amplitude(s) for diagram number 7173
    FFV1_0( w_fp[355], w_fp[481], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7174 OF 15495 ***
    // Wavefunction(s) for diagram number 7174
    FFV1_1( w_fp[449], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[438] );
    // Amplitude(s) for diagram number 7174
    FFV1_0( w_fp[216], w_fp[438], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7175 OF 15495 ***
    // Wavefunction(s) for diagram number 7175
    // (none)
    // Amplitude(s) for diagram number 7175
    FFV1_0( w_fp[216], w_fp[481], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7176 OF 15495 ***
    // Wavefunction(s) for diagram number 7176
    // (none)
    // Amplitude(s) for diagram number 7176
    FFV1_0( w_fp[199], w_fp[438], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7177 OF 15495 ***
    // Wavefunction(s) for diagram number 7177
    // (none)
    // Amplitude(s) for diagram number 7177
    FFV1_0( w_fp[199], w_fp[590], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7178 OF 15495 ***
    // Wavefunction(s) for diagram number 7178
    // (none)
    // Amplitude(s) for diagram number 7178
    FFV1_0( w_fp[355], w_fp[594], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];

    // *** DIAGRAM 7179 OF 15495 ***
    // Wavefunction(s) for diagram number 7179
    // (none)
    // Amplitude(s) for diagram number 7179
    FFV1_0( w_fp[355], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7180 OF 15495 ***
    // Wavefunction(s) for diagram number 7180
    // (none)
    // Amplitude(s) for diagram number 7180
    VVVV1_0( w_fp[474], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] += amp_sv[0];
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

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
#endif
#endif
  }

  //--------------------------------------------------------------------------


#ifndef MGONGPU_RDC_DIAGRAMS
  __global__ void
#else
  __device__ void
#endif
  diagramgroup719( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 528 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 594 );
#endif
#endif

    // *** DIAGRAM 7181 OF 15495 ***
    // Wavefunction(s) for diagram number 7181
    // (none)
    // Amplitude(s) for diagram number 7181
    VVV1_0( w_fp[239], w_fp[7], w_fp[545], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[640] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 7182 OF 15495 ***
    // Wavefunction(s) for diagram number 7182
    // (none)
    // Amplitude(s) for diagram number 7182
    VVV1_0( w_fp[1], w_fp[239], w_fp[536], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 7183 OF 15495 ***
    // Wavefunction(s) for diagram number 7183
    // (none)
    // Amplitude(s) for diagram number 7183
    FFV1_0( w_fp[199], w_fp[2], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7184 OF 15495 ***
    // Wavefunction(s) for diagram number 7184
    // (none)
    // Amplitude(s) for diagram number 7184
    FFV1_0( w_fp[199], w_fp[594], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];

    // *** DIAGRAM 7185 OF 15495 ***
    // Wavefunction(s) for diagram number 7185
    // (none)
    // Amplitude(s) for diagram number 7185
    FFV1_0( w_fp[355], w_fp[571], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 7186 OF 15495 ***
    // Wavefunction(s) for diagram number 7186
    // (none)
    // Amplitude(s) for diagram number 7186
    FFV1_0( w_fp[355], w_fp[2], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7187 OF 15495 ***
    // Wavefunction(s) for diagram number 7187
    // (none)
    // Amplitude(s) for diagram number 7187
    VVVV1_0( w_fp[547], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    VVVV3_0( w_fp[547], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    VVVV4_0( w_fp[547], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 7188 OF 15495 ***
    // Wavefunction(s) for diagram number 7188
    // (none)
    // Amplitude(s) for diagram number 7188
    VVV1_0( w_fp[239], w_fp[5], w_fp[528], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7189 OF 15495 ***
    // Wavefunction(s) for diagram number 7189
    // (none)
    // Amplitude(s) for diagram number 7189
    VVV1_0( w_fp[1], w_fp[239], w_fp[519], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];

    // *** DIAGRAM 7190 OF 15495 ***
    // Wavefunction(s) for diagram number 7190
    // (none)
    // Amplitude(s) for diagram number 7190
    FFV1_0( w_fp[216], w_fp[2], w_fp[528], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];

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
  diagramgroup720( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 589 );
#endif
#endif

    // *** DIAGRAM 7191 OF 15495 ***
    // Wavefunction(s) for diagram number 7191
    // (none)
    // Amplitude(s) for diagram number 7191
    FFV1_0( w_fp[216], w_fp[571], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];

    // *** DIAGRAM 7192 OF 15495 ***
    // Wavefunction(s) for diagram number 7192
    // (none)
    // Amplitude(s) for diagram number 7192
    VVV1_0( w_fp[525], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    VVV1_0( w_fp[526], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[640] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    VVV1_0( w_fp[45], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[682] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];

    // *** DIAGRAM 7193 OF 15495 ***
    // Wavefunction(s) for diagram number 7193
    // (none)
    // Amplitude(s) for diagram number 7193
    FFV1_0( w_fp[199], w_fp[2], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[85] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[45], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7194 OF 15495 ***
    // Wavefunction(s) for diagram number 7194
    // (none)
    // Amplitude(s) for diagram number 7194
    VVV1_0( w_fp[583], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += amp_sv[0];
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[466] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    VVV1_0( w_fp[584], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    VVV1_0( w_fp[561], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];

    // *** DIAGRAM 7195 OF 15495 ***
    // Wavefunction(s) for diagram number 7195
    // (none)
    // Amplitude(s) for diagram number 7195
    FFV1_0( w_fp[216], w_fp[2], w_fp[583], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[584], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7196 OF 15495 ***
    // Wavefunction(s) for diagram number 7196
    // (none)
    // Amplitude(s) for diagram number 7196
    FFV1_0( w_fp[355], w_fp[2], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7197 OF 15495 ***
    // Wavefunction(s) for diagram number 7197
    // (none)
    // Amplitude(s) for diagram number 7197
    VVV1_0( w_fp[454], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[135] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[232] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    VVV1_0( w_fp[575], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[94] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    VVV1_0( w_fp[589], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[88] += amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[232] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 7198 OF 15495 ***
    // Wavefunction(s) for diagram number 7198
    // (none)
    // Amplitude(s) for diagram number 7198
    VVV1_0( w_fp[450], w_fp[239], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[135] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[466] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    jamp_sv[682] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];

    // *** DIAGRAM 7199 OF 15495 ***
    // Wavefunction(s) for diagram number 7199
    // (none)
    // Amplitude(s) for diagram number 7199
    FFV1_0( w_fp[196], w_fp[122], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 7200 OF 15495 ***
    // Wavefunction(s) for diagram number 7200
    // (none)
    // Amplitude(s) for diagram number 7200
    FFV1_0( w_fp[123], w_fp[2], w_fp[450], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[135] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];

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
