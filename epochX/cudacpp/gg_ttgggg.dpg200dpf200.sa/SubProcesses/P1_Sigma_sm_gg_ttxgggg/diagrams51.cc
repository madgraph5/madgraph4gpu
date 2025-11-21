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
  diagramgroup51( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 112 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 184 );
    retrieveWf( wfs, w_cx, nevt, 186 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 220 );
    retrieveWf( wfs, w_cx, nevt, 224 );
    retrieveWf( wfs, w_cx, nevt, 226 );
    retrieveWf( wfs, w_cx, nevt, 227 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 231 );
    retrieveWf( wfs, w_cx, nevt, 232 );
    retrieveWf( wfs, w_cx, nevt, 235 );
    retrieveWf( wfs, w_cx, nevt, 258 );
    retrieveWf( wfs, w_cx, nevt, 260 );
    retrieveWf( wfs, w_cx, nevt, 263 );
    retrieveWf( wfs, w_cx, nevt, 265 );
    retrieveWf( wfs, w_cx, nevt, 266 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 269 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 271 );
    retrieveWf( wfs, w_cx, nevt, 272 );
    retrieveWf( wfs, w_cx, nevt, 273 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 275 );
    retrieveWf( wfs, w_cx, nevt, 281 );
    retrieveWf( wfs, w_cx, nevt, 283 );
    retrieveWf( wfs, w_cx, nevt, 286 );
    retrieveWf( wfs, w_cx, nevt, 287 );
    retrieveWf( wfs, w_cx, nevt, 288 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 483 );
    retrieveWf( wfs, w_cx, nevt, 491 );
    retrieveWf( wfs, w_cx, nevt, 493 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 542 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 680 );
    retrieveWf( wfs, w_cx, nevt, 690 );
    retrieveWf( wfs, w_cx, nevt, 691 );
    retrieveWf( wfs, w_cx, nevt, 693 );
    retrieveWf( wfs, w_cx, nevt, 694 );
    retrieveWf( wfs, w_cx, nevt, 700 );
    retrieveWf( wfs, w_cx, nevt, 703 );
    retrieveWf( wfs, w_cx, nevt, 705 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 707 );
    retrieveWf( wfs, w_cx, nevt, 708 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 711 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
    retrieveWf( wfs, w_cx, nevt, 724 );
    retrieveWf( wfs, w_cx, nevt, 725 );
    retrieveWf( wfs, w_cx, nevt, 726 );
    retrieveWf( wfs, w_cx, nevt, 727 );
    retrieveWf( wfs, w_cx, nevt, 728 );
    retrieveWf( wfs, w_cx, nevt, 729 );
    retrieveWf( wfs, w_cx, nevt, 733 );
    retrieveWf( wfs, w_cx, nevt, 734 );
    retrieveWf( wfs, w_cx, nevt, 735 );
    retrieveWf( wfs, w_cx, nevt, 741 );
    retrieveWf( wfs, w_cx, nevt, 743 );
#endif
#endif

    // *** DIAGRAM 10001 OF 15495 ***
    // Wavefunction(s) for diagram number 10001
    // (none)
    // Amplitude(s) for diagram number 10001
    VVV1_0( w_fp[0], w_fp[260], w_fp[220], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10002 OF 15495 ***
    // Wavefunction(s) for diagram number 10002
    // (none)
    // Amplitude(s) for diagram number 10002
    FFV1_0( w_fp[179], w_fp[197], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[197], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[197], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[552] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10003 OF 15495 ***
    // Wavefunction(s) for diagram number 10003
    // (none)
    // Amplitude(s) for diagram number 10003
    VVV1_0( w_fp[662], w_fp[214], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 10004 OF 15495 ***
    // Wavefunction(s) for diagram number 10004
    // (none)
    // Amplitude(s) for diagram number 10004
    FFV1_0( w_fp[3], w_fp[224], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[573] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[575] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[597] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[599] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10005 OF 15495 ***
    // Wavefunction(s) for diagram number 10005
    // (none)
    // Amplitude(s) for diagram number 10005
    FFV1_0( w_fp[186], w_fp[197], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10006 OF 15495 ***
    // Wavefunction(s) for diagram number 10006
    // (none)
    // Amplitude(s) for diagram number 10006
    FFV1_0( w_fp[275], w_fp[540], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];

    // *** DIAGRAM 10007 OF 15495 ***
    // Wavefunction(s) for diagram number 10007
    // (none)
    // Amplitude(s) for diagram number 10007
    FFV1_0( w_fp[3], w_fp[540], w_fp[287], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10008 OF 15495 ***
    // Wavefunction(s) for diagram number 10008
    // (none)
    // Amplitude(s) for diagram number 10008
    FFV1_0( w_fp[186], w_fp[540], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];

    // *** DIAGRAM 10009 OF 15495 ***
    // Wavefunction(s) for diagram number 10009
    // (none)
    // Amplitude(s) for diagram number 10009
    FFV1_0( w_fp[3], w_fp[491], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[513] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[515] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[537] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[539] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10010 OF 15495 ***
    // Wavefunction(s) for diagram number 10010
    // (none)
    // Amplitude(s) for diagram number 10010
    FFV1_0( w_fp[275], w_fp[197], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[497] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[571] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[595] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10011 OF 15495 ***
    // Wavefunction(s) for diagram number 10011
    // (none)
    // Amplitude(s) for diagram number 10011
    VVV1_0( w_fp[522], w_fp[258], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[497] += amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[515] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[537] += amp_sv[0];
    jamp_sv[539] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];

    // *** DIAGRAM 10012 OF 15495 ***
    // Wavefunction(s) for diagram number 10012
    // (none)
    // Amplitude(s) for diagram number 10012
    FFV1_0( w_fp[186], w_fp[491], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[510] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 10013 OF 15495 ***
    // Wavefunction(s) for diagram number 10013
    // (none)
    // Amplitude(s) for diagram number 10013
    FFV1_0( w_fp[275], w_fp[224], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += amp_sv[0];
    jamp_sv[571] -= amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[595] += amp_sv[0];

    // *** DIAGRAM 10014 OF 15495 ***
    // Wavefunction(s) for diagram number 10014
    // (none)
    // Amplitude(s) for diagram number 10014
    VVV1_0( w_fp[0], w_fp[287], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];

    // *** DIAGRAM 10015 OF 15495 ***
    // Wavefunction(s) for diagram number 10015
    // (none)
    // Amplitude(s) for diagram number 10015
    FFV1_0( w_fp[3], w_fp[197], w_fp[733], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[573] -= amp_sv[0];
    jamp_sv[575] += amp_sv[0];
    jamp_sv[597] += amp_sv[0];
    jamp_sv[599] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[734], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += amp_sv[0];
    jamp_sv[497] -= amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[515] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[539] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[735], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[571] += amp_sv[0];
    jamp_sv[573] += amp_sv[0];
    jamp_sv[575] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[595] -= amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[599] += amp_sv[0];

    // *** DIAGRAM 10016 OF 15495 ***
    // Wavefunction(s) for diagram number 10016
    // (none)
    // Amplitude(s) for diagram number 10016
    VVVV1_0( w_fp[662], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[228], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 10017 OF 15495 ***
    // Wavefunction(s) for diagram number 10017
    // (none)
    // Amplitude(s) for diagram number 10017
    VVV1_0( w_fp[228], w_fp[6], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10018 OF 15495 ***
    // Wavefunction(s) for diagram number 10018
    // (none)
    // Amplitude(s) for diagram number 10018
    VVV1_0( w_fp[228], w_fp[5], w_fp[670], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 10019 OF 15495 ***
    // Wavefunction(s) for diagram number 10019
    // (none)
    // Amplitude(s) for diagram number 10019
    FFV1_0( w_fp[741], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];

    // *** DIAGRAM 10020 OF 15495 ***
    // Wavefunction(s) for diagram number 10020
    // (none)
    // Amplitude(s) for diagram number 10020
    FFV1_0( w_fp[3], w_fp[226], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10021 OF 15495 ***
    // Wavefunction(s) for diagram number 10021
    // (none)
    // Amplitude(s) for diagram number 10021
    FFV1_0( w_fp[741], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10022 OF 15495 ***
    // Wavefunction(s) for diagram number 10022
    // (none)
    // Amplitude(s) for diagram number 10022
    FFV1_0( w_fp[3], w_fp[227], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10023 OF 15495 ***
    // Wavefunction(s) for diagram number 10023
    // (none)
    // Amplitude(s) for diagram number 10023
    FFV1_0( w_fp[275], w_fp[690], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10024 OF 15495 ***
    // Wavefunction(s) for diagram number 10024
    // (none)
    // Amplitude(s) for diagram number 10024
    FFV1_0( w_fp[275], w_fp[691], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10025 OF 15495 ***
    // Wavefunction(s) for diagram number 10025
    FFV1P0_3( w_fp[3], w_fp[542], COUPs[1], 1.0, depCoup, 0., 0., w_fp[491] );
    // Amplitude(s) for diagram number 10025
    VVV1_0( w_fp[260], w_fp[6], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10026 OF 15495 ***
    // Wavefunction(s) for diagram number 10026
    // (none)
    // Amplitude(s) for diagram number 10026
    FFV1_0( w_fp[3], w_fp[691], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];

    // *** DIAGRAM 10027 OF 15495 ***
    // Wavefunction(s) for diagram number 10027
    // (none)
    // Amplitude(s) for diagram number 10027
    VVV1_0( w_fp[263], w_fp[5], w_fp[491], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10028 OF 15495 ***
    // Wavefunction(s) for diagram number 10028
    // (none)
    // Amplitude(s) for diagram number 10028
    FFV1_0( w_fp[3], w_fp[690], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];

    // *** DIAGRAM 10029 OF 15495 ***
    // Wavefunction(s) for diagram number 10029
    // (none)
    // Amplitude(s) for diagram number 10029
    FFV1_0( w_fp[3], w_fp[542], w_fp[266], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[267], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[542], w_fp[268], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10030 OF 15495 ***
    // Wavefunction(s) for diagram number 10030
    // (none)
    // Amplitude(s) for diagram number 10030
    FFV1_0( w_fp[743], w_fp[226], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10031 OF 15495 ***
    // Wavefunction(s) for diagram number 10031
    // (none)
    // Amplitude(s) for diagram number 10031
    FFV1_0( w_fp[275], w_fp[693], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10032 OF 15495 ***
    // Wavefunction(s) for diagram number 10032
    // (none)
    // Amplitude(s) for diagram number 10032
    FFV1_0( w_fp[743], w_fp[227], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10033 OF 15495 ***
    // Wavefunction(s) for diagram number 10033
    // (none)
    // Amplitude(s) for diagram number 10033
    FFV1_0( w_fp[275], w_fp[694], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10034 OF 15495 ***
    // Wavefunction(s) for diagram number 10034
    // (none)
    // Amplitude(s) for diagram number 10034
    VVVV1_0( w_fp[0], w_fp[260], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[260], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[260], w_fp[228], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 10035 OF 15495 ***
    // Wavefunction(s) for diagram number 10035
    // (none)
    // Amplitude(s) for diagram number 10035
    VVV1_0( w_fp[228], w_fp[6], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10036 OF 15495 ***
    // Wavefunction(s) for diagram number 10036
    VVV1P0_1( w_fp[0], w_fp[228], COUPs[0], 1.0, depCoup, 0., 0., w_fp[746] );
    // Amplitude(s) for diagram number 10036
    VVV1_0( w_fp[260], w_fp[6], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10037 OF 15495 ***
    // Wavefunction(s) for diagram number 10037
    // (none)
    // Amplitude(s) for diagram number 10037
    FFV1_0( w_fp[3], w_fp[227], w_fp[700], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10038 OF 15495 ***
    // Wavefunction(s) for diagram number 10038
    // (none)
    // Amplitude(s) for diagram number 10038
    FFV1_0( w_fp[3], w_fp[694], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[701] += amp_sv[0];

    // *** DIAGRAM 10039 OF 15495 ***
    // Wavefunction(s) for diagram number 10039
    // (none)
    // Amplitude(s) for diagram number 10039
    VVVV1_0( w_fp[0], w_fp[263], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[263], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[263], w_fp[228], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 10040 OF 15495 ***
    // Wavefunction(s) for diagram number 10040
    // (none)
    // Amplitude(s) for diagram number 10040
    VVV1_0( w_fp[228], w_fp[5], w_fp[703], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 10041 OF 15495 ***
    // Wavefunction(s) for diagram number 10041
    // (none)
    // Amplitude(s) for diagram number 10041
    VVV1_0( w_fp[263], w_fp[5], w_fp[746], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 10042 OF 15495 ***
    // Wavefunction(s) for diagram number 10042
    // (none)
    // Amplitude(s) for diagram number 10042
    FFV1_0( w_fp[3], w_fp[226], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10043 OF 15495 ***
    // Wavefunction(s) for diagram number 10043
    // (none)
    // Amplitude(s) for diagram number 10043
    FFV1_0( w_fp[3], w_fp[693], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];

    // *** DIAGRAM 10044 OF 15495 ***
    // Wavefunction(s) for diagram number 10044
    // (none)
    // Amplitude(s) for diagram number 10044
    VVV1_0( w_fp[706], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];
    VVV1_0( w_fp[707], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[701] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    VVV1_0( w_fp[708], w_fp[228], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 10045 OF 15495 ***
    // Wavefunction(s) for diagram number 10045
    // (none)
    // Amplitude(s) for diagram number 10045
    FFV1_0( w_fp[3], w_fp[227], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[701] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[705] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[227], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10046 OF 15495 ***
    // Wavefunction(s) for diagram number 10046
    // (none)
    // Amplitude(s) for diagram number 10046
    VVV1_0( w_fp[709], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0( w_fp[710], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    VVV1_0( w_fp[711], w_fp[228], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10047 OF 15495 ***
    // Wavefunction(s) for diagram number 10047
    // (none)
    // Amplitude(s) for diagram number 10047
    FFV1_0( w_fp[3], w_fp[226], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[681] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[226], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10048 OF 15495 ***
    // Wavefunction(s) for diagram number 10048
    // (none)
    // Amplitude(s) for diagram number 10048
    VVV1_0( w_fp[0], w_fp[266], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[267], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[268], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[687] -= amp_sv[0];
    jamp_sv[705] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10049 OF 15495 ***
    // Wavefunction(s) for diagram number 10049
    // (none)
    // Amplitude(s) for diagram number 10049
    VVV1_0( w_fp[662], w_fp[231], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10050 OF 15495 ***
    // Wavefunction(s) for diagram number 10050
    // (none)
    // Amplitude(s) for diagram number 10050
    FFV1_0( w_fp[168], w_fp[227], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10051 OF 15495 ***
    // Wavefunction(s) for diagram number 10051
    // (none)
    // Amplitude(s) for diagram number 10051
    FFV1_0( w_fp[170], w_fp[215], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 10052 OF 15495 ***
    // Wavefunction(s) for diagram number 10052
    // (none)
    // Amplitude(s) for diagram number 10052
    FFV1_0( w_fp[281], w_fp[542], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10053 OF 15495 ***
    // Wavefunction(s) for diagram number 10053
    // (none)
    // Amplitude(s) for diagram number 10053
    FFV1_0( w_fp[168], w_fp[542], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 10054 OF 15495 ***
    // Wavefunction(s) for diagram number 10054
    // (none)
    // Amplitude(s) for diagram number 10054
    FFV1_0( w_fp[170], w_fp[542], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10055 OF 15495 ***
    // Wavefunction(s) for diagram number 10055
    // (none)
    // Amplitude(s) for diagram number 10055
    FFV1_0( w_fp[527], w_fp[493], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10056 OF 15495 ***
    // Wavefunction(s) for diagram number 10056
    // (none)
    // Amplitude(s) for diagram number 10056
    FFV1_0( w_fp[527], w_fp[215], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10057 OF 15495 ***
    // Wavefunction(s) for diagram number 10057
    // (none)
    // Amplitude(s) for diagram number 10057
    FFV1_0( w_fp[527], w_fp[227], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10058 OF 15495 ***
    // Wavefunction(s) for diagram number 10058
    // (none)
    // Amplitude(s) for diagram number 10058
    FFV1_0( w_fp[170], w_fp[493], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10059 OF 15495 ***
    // Wavefunction(s) for diagram number 10059
    // (none)
    // Amplitude(s) for diagram number 10059
    FFV1_0( w_fp[281], w_fp[227], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10060 OF 15495 ***
    // Wavefunction(s) for diagram number 10060
    // (none)
    // Amplitude(s) for diagram number 10060
    VVV1_0( w_fp[0], w_fp[263], w_fp[231], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10061 OF 15495 ***
    // Wavefunction(s) for diagram number 10061
    // (none)
    // Amplitude(s) for diagram number 10061
    FFV1_0( w_fp[168], w_fp[215], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[215], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10062 OF 15495 ***
    // Wavefunction(s) for diagram number 10062
    // (none)
    // Amplitude(s) for diagram number 10062
    VVV1_0( w_fp[662], w_fp[232], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10063 OF 15495 ***
    // Wavefunction(s) for diagram number 10063
    // (none)
    // Amplitude(s) for diagram number 10063
    FFV1_0( w_fp[174], w_fp[226], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10064 OF 15495 ***
    // Wavefunction(s) for diagram number 10064
    // (none)
    // Amplitude(s) for diagram number 10064
    FFV1_0( w_fp[175], w_fp[215], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 10065 OF 15495 ***
    // Wavefunction(s) for diagram number 10065
    // (none)
    // Amplitude(s) for diagram number 10065
    FFV1_0( w_fp[283], w_fp[542], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10066 OF 15495 ***
    // Wavefunction(s) for diagram number 10066
    // (none)
    // Amplitude(s) for diagram number 10066
    FFV1_0( w_fp[174], w_fp[542], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];

    // *** DIAGRAM 10067 OF 15495 ***
    // Wavefunction(s) for diagram number 10067
    // (none)
    // Amplitude(s) for diagram number 10067
    FFV1_0( w_fp[175], w_fp[542], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10068 OF 15495 ***
    // Wavefunction(s) for diagram number 10068
    // (none)
    // Amplitude(s) for diagram number 10068
    FFV1_0( w_fp[594], w_fp[493], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10069 OF 15495 ***
    // Wavefunction(s) for diagram number 10069
    // (none)
    // Amplitude(s) for diagram number 10069
    FFV1_0( w_fp[594], w_fp[215], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10070 OF 15495 ***
    // Wavefunction(s) for diagram number 10070
    // (none)
    // Amplitude(s) for diagram number 10070
    FFV1_0( w_fp[594], w_fp[226], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10071 OF 15495 ***
    // Wavefunction(s) for diagram number 10071
    // (none)
    // Amplitude(s) for diagram number 10071
    FFV1_0( w_fp[175], w_fp[493], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10072 OF 15495 ***
    // Wavefunction(s) for diagram number 10072
    // (none)
    // Amplitude(s) for diagram number 10072
    FFV1_0( w_fp[283], w_fp[226], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10073 OF 15495 ***
    // Wavefunction(s) for diagram number 10073
    // (none)
    // Amplitude(s) for diagram number 10073
    VVV1_0( w_fp[0], w_fp[260], w_fp[232], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10074 OF 15495 ***
    // Wavefunction(s) for diagram number 10074
    // (none)
    // Amplitude(s) for diagram number 10074
    FFV1_0( w_fp[174], w_fp[215], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10075 OF 15495 ***
    // Wavefunction(s) for diagram number 10075
    // (none)
    // Amplitude(s) for diagram number 10075
    VVV1_0( w_fp[662], w_fp[228], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10076 OF 15495 ***
    // Wavefunction(s) for diagram number 10076
    // (none)
    // Amplitude(s) for diagram number 10076
    FFV1_0( w_fp[3], w_fp[235], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[693] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[695] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[717] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[719] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10077 OF 15495 ***
    // Wavefunction(s) for diagram number 10077
    // (none)
    // Amplitude(s) for diagram number 10077
    FFV1_0( w_fp[184], w_fp[215], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10078 OF 15495 ***
    // Wavefunction(s) for diagram number 10078
    // (none)
    // Amplitude(s) for diagram number 10078
    FFV1_0( w_fp[275], w_fp[542], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];

    // *** DIAGRAM 10079 OF 15495 ***
    // Wavefunction(s) for diagram number 10079
    // (none)
    // Amplitude(s) for diagram number 10079
    FFV1_0( w_fp[3], w_fp[542], w_fp[286], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10080 OF 15495 ***
    // Wavefunction(s) for diagram number 10080
    // (none)
    // Amplitude(s) for diagram number 10080
    FFV1_0( w_fp[184], w_fp[542], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];

    // *** DIAGRAM 10081 OF 15495 ***
    // Wavefunction(s) for diagram number 10081
    // (none)
    // Amplitude(s) for diagram number 10081
    FFV1_0( w_fp[3], w_fp[493], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[633] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[657] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[659] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10082 OF 15495 ***
    // Wavefunction(s) for diagram number 10082
    // (none)
    // Amplitude(s) for diagram number 10082
    FFV1_0( w_fp[275], w_fp[215], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[623] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[690] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[691] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[715] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10083 OF 15495 ***
    // Wavefunction(s) for diagram number 10083
    // (none)
    // Amplitude(s) for diagram number 10083
    VVV1_0( w_fp[559], w_fp[258], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    jamp_sv[623] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[657] += amp_sv[0];
    jamp_sv[659] -= amp_sv[0];
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 10084 OF 15495 ***
    // Wavefunction(s) for diagram number 10084
    // (none)
    // Amplitude(s) for diagram number 10084
    FFV1_0( w_fp[184], w_fp[493], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 10085 OF 15495 ***
    // Wavefunction(s) for diagram number 10085
    // (none)
    // Amplitude(s) for diagram number 10085
    FFV1_0( w_fp[275], w_fp[235], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[690] += amp_sv[0];
    jamp_sv[691] -= amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    jamp_sv[715] += amp_sv[0];

    // *** DIAGRAM 10086 OF 15495 ***
    // Wavefunction(s) for diagram number 10086
    // (none)
    // Amplitude(s) for diagram number 10086
    VVV1_0( w_fp[0], w_fp[286], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];

    // *** DIAGRAM 10087 OF 15495 ***
    // Wavefunction(s) for diagram number 10087
    // (none)
    // Amplitude(s) for diagram number 10087
    FFV1_0( w_fp[3], w_fp[215], w_fp[727], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[693] -= amp_sv[0];
    jamp_sv[695] += amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    jamp_sv[719] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[728], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[616] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    jamp_sv[623] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[633] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[659] += amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[729], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[600] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[690] -= amp_sv[0];
    jamp_sv[691] += amp_sv[0];
    jamp_sv[693] += amp_sv[0];
    jamp_sv[695] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    jamp_sv[715] -= amp_sv[0];
    jamp_sv[717] -= amp_sv[0];
    jamp_sv[719] += amp_sv[0];

    // *** DIAGRAM 10088 OF 15495 ***
    // Wavefunction(s) for diagram number 10088
    // (none)
    // Amplitude(s) for diagram number 10088
    VVVV1_0( w_fp[662], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[154], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 10089 OF 15495 ***
    // Wavefunction(s) for diagram number 10089
    // (none)
    // Amplitude(s) for diagram number 10089
    VVV1_0( w_fp[154], w_fp[7], w_fp[670], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10090 OF 15495 ***
    // Wavefunction(s) for diagram number 10090
    // (none)
    // Amplitude(s) for diagram number 10090
    VVV1_0( w_fp[154], w_fp[6], w_fp[680], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 10091 OF 15495 ***
    // Wavefunction(s) for diagram number 10091
    FFV1_1( w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[493] );
    // Amplitude(s) for diagram number 10091
    FFV1_0( w_fp[170], w_fp[493], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];

    // *** DIAGRAM 10092 OF 15495 ***
    // Wavefunction(s) for diagram number 10092
    // (none)
    // Amplitude(s) for diagram number 10092
    FFV1_0( w_fp[170], w_fp[2], w_fp[680], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10093 OF 15495 ***
    // Wavefunction(s) for diagram number 10093
    // (none)
    // Amplitude(s) for diagram number 10093
    FFV1_0( w_fp[171], w_fp[493], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];

    // *** DIAGRAM 10094 OF 15495 ***
    // Wavefunction(s) for diagram number 10094
    // (none)
    // Amplitude(s) for diagram number 10094
    FFV1_0( w_fp[171], w_fp[2], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10095 OF 15495 ***
    // Wavefunction(s) for diagram number 10095
    // (none)
    // Amplitude(s) for diagram number 10095
    FFV1_0( w_fp[538], w_fp[483], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10096 OF 15495 ***
    // Wavefunction(s) for diagram number 10096
    // (none)
    // Amplitude(s) for diagram number 10096
    FFV1_0( w_fp[543], w_fp[483], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10097 OF 15495 ***
    // Wavefunction(s) for diagram number 10097
    FFV1P0_3( w_fp[527], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[747] );
    // Amplitude(s) for diagram number 10097
    VVV1_0( w_fp[263], w_fp[7], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10098 OF 15495 ***
    // Wavefunction(s) for diagram number 10098
    // (none)
    // Amplitude(s) for diagram number 10098
    FFV1_0( w_fp[543], w_fp[2], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];

    // *** DIAGRAM 10099 OF 15495 ***
    // Wavefunction(s) for diagram number 10099
    // (none)
    // Amplitude(s) for diagram number 10099
    VVV1_0( w_fp[265], w_fp[6], w_fp[747], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10100 OF 15495 ***
    // Wavefunction(s) for diagram number 10100
    // (none)
    // Amplitude(s) for diagram number 10100
    FFV1_0( w_fp[538], w_fp[2], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[658] += amp_sv[0];

    // *** DIAGRAM 10101 OF 15495 ***
    // Wavefunction(s) for diagram number 10101
    // (none)
    // Amplitude(s) for diagram number 10101
    FFV1_0( w_fp[527], w_fp[2], w_fp[272], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[273], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[527], w_fp[2], w_fp[274], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10102 OF 15495 ***
    // Wavefunction(s) for diagram number 10102
    FFV1_1( w_fp[483], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[748] );
    // Amplitude(s) for diagram number 10102
    FFV1_0( w_fp[170], w_fp[748], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10103 OF 15495 ***
    // Wavefunction(s) for diagram number 10103
    // (none)
    // Amplitude(s) for diagram number 10103
    FFV1_0( w_fp[535], w_fp[483], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10104 OF 15495 ***
    // Wavefunction(s) for diagram number 10104
    // (none)
    // Amplitude(s) for diagram number 10104
    FFV1_0( w_fp[171], w_fp[748], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10105 OF 15495 ***
    // Wavefunction(s) for diagram number 10105
    // (none)
    // Amplitude(s) for diagram number 10105
    FFV1_0( w_fp[450], w_fp[483], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10106 OF 15495 ***
    // Wavefunction(s) for diagram number 10106
    // (none)
    // Amplitude(s) for diagram number 10106
    VVVV1_0( w_fp[0], w_fp[263], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[263], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[263], w_fp[154], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 10107 OF 15495 ***
    // Wavefunction(s) for diagram number 10107
    // (none)
    // Amplitude(s) for diagram number 10107
    VVV1_0( w_fp[154], w_fp[7], w_fp[703], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10108 OF 15495 ***
    // Wavefunction(s) for diagram number 10108
    VVV1P0_1( w_fp[0], w_fp[154], COUPs[0], 1.0, depCoup, 0., 0., w_fp[749] );
    // Amplitude(s) for diagram number 10108
    VVV1_0( w_fp[263], w_fp[7], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10109 OF 15495 ***
    // Wavefunction(s) for diagram number 10109
    // (none)
    // Amplitude(s) for diagram number 10109
    FFV1_0( w_fp[171], w_fp[2], w_fp[703], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10110 OF 15495 ***
    // Wavefunction(s) for diagram number 10110
    // (none)
    // Amplitude(s) for diagram number 10110
    FFV1_0( w_fp[450], w_fp[2], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];

    // *** DIAGRAM 10111 OF 15495 ***
    // Wavefunction(s) for diagram number 10111
    // (none)
    // Amplitude(s) for diagram number 10111
    VVVV1_0( w_fp[0], w_fp[265], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[265], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[265], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 10112 OF 15495 ***
    // Wavefunction(s) for diagram number 10112
    // (none)
    // Amplitude(s) for diagram number 10112
    VVV1_0( w_fp[154], w_fp[6], w_fp[705], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];

    // *** DIAGRAM 10113 OF 15495 ***
    // Wavefunction(s) for diagram number 10113
    // (none)
    // Amplitude(s) for diagram number 10113
    VVV1_0( w_fp[265], w_fp[6], w_fp[749], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];

    // *** DIAGRAM 10114 OF 15495 ***
    // Wavefunction(s) for diagram number 10114
    // (none)
    // Amplitude(s) for diagram number 10114
    FFV1_0( w_fp[170], w_fp[2], w_fp[705], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10115 OF 15495 ***
    // Wavefunction(s) for diagram number 10115
    // (none)
    // Amplitude(s) for diagram number 10115
    FFV1_0( w_fp[535], w_fp[2], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[163] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 10116 OF 15495 ***
    // Wavefunction(s) for diagram number 10116
    // (none)
    // Amplitude(s) for diagram number 10116
    VVV1_0( w_fp[709], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    VVV1_0( w_fp[710], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    VVV1_0( w_fp[711], w_fp[154], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 10117 OF 15495 ***
    // Wavefunction(s) for diagram number 10117
    // (none)
    // Amplitude(s) for diagram number 10117
    FFV1_0( w_fp[171], w_fp[2], w_fp[709], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[710], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[73] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[157] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[171], w_fp[2], w_fp[711], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10118 OF 15495 ***
    // Wavefunction(s) for diagram number 10118
    // (none)
    // Amplitude(s) for diagram number 10118
    VVV1_0( w_fp[712], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[481] -= amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[655] -= amp_sv[0];
    VVV1_0( w_fp[713], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    VVV1_0( w_fp[714], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[481] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[655] += amp_sv[0];

    // *** DIAGRAM 10119 OF 15495 ***
    // Wavefunction(s) for diagram number 10119
    // (none)
    // Amplitude(s) for diagram number 10119
    FFV1_0( w_fp[170], w_fp[2], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[713], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[714], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[655] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10120 OF 15495 ***
    // Wavefunction(s) for diagram number 10120
    // (none)
    // Amplitude(s) for diagram number 10120
    VVV1_0( w_fp[0], w_fp[272], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[273], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[5] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[274], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[73] -= amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[103] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[658] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10121 OF 15495 ***
    // Wavefunction(s) for diagram number 10121
    // (none)
    // Amplitude(s) for diagram number 10121
    VVV1_0( w_fp[662], w_fp[154], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10122 OF 15495 ***
    // Wavefunction(s) for diagram number 10122
    // (none)
    // Amplitude(s) for diagram number 10122
    FFV1_0( w_fp[168], w_fp[128], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10123 OF 15495 ***
    // Wavefunction(s) for diagram number 10123
    // (none)
    // Amplitude(s) for diagram number 10123
    FFV1_0( w_fp[112], w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[5] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10124 OF 15495 ***
    // Wavefunction(s) for diagram number 10124
    // (none)
    // Amplitude(s) for diagram number 10124
    FFV1_0( w_fp[527], w_fp[483], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];

    // *** DIAGRAM 10125 OF 15495 ***
    // Wavefunction(s) for diagram number 10125
    // (none)
    // Amplitude(s) for diagram number 10125
    FFV1_0( w_fp[527], w_fp[2], w_fp[288], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10126 OF 15495 ***
    // Wavefunction(s) for diagram number 10126
    // (none)
    // Amplitude(s) for diagram number 10126
    FFV1_0( w_fp[527], w_fp[128], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10127 OF 15495 ***
    // Wavefunction(s) for diagram number 10127
    // (none)
    // Amplitude(s) for diagram number 10127
    FFV1_0( w_fp[168], w_fp[483], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[267] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[280] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[286] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10128 OF 15495 ***
    // Wavefunction(s) for diagram number 10128
    // (none)
    // Amplitude(s) for diagram number 10128
    FFV1_0( w_fp[281], w_fp[2], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[92] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[696] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10129 OF 15495 ***
    // Wavefunction(s) for diagram number 10129
    // (none)
    // Amplitude(s) for diagram number 10129
    VVV1_0( w_fp[516], w_fp[258], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[286] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];

    // *** DIAGRAM 10130 OF 15495 ***
    // Wavefunction(s) for diagram number 10130
    // (none)
    // Amplitude(s) for diagram number 10130
    FFV1_0( w_fp[112], w_fp[483], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[267] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];

    // *** DIAGRAM 10131 OF 15495 ***
    // Wavefunction(s) for diagram number 10131
    // (none)
    // Amplitude(s) for diagram number 10131
    FFV1_0( w_fp[281], w_fp[128], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[696] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];

    // *** DIAGRAM 10132 OF 15495 ***
    // Wavefunction(s) for diagram number 10132
    // (none)
    // Amplitude(s) for diagram number 10132
    VVV1_0( w_fp[0], w_fp[288], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];

    // *** DIAGRAM 10133 OF 15495 ***
    // Wavefunction(s) for diagram number 10133
    // (none)
    // Amplitude(s) for diagram number 10133
    FFV1_0( w_fp[168], w_fp[2], w_fp[726], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] -= amp_sv[0];
    jamp_sv[5] += amp_sv[0];
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[725], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[90] += amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[286] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[724], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[3] += amp_sv[0];
    jamp_sv[5] -= amp_sv[0];
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[267] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[590] -= amp_sv[0];
    jamp_sv[696] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[710] += amp_sv[0];

    // *** DIAGRAM 10134 OF 15495 ***
    // Wavefunction(s) for diagram number 10134
    // (none)
    // Amplitude(s) for diagram number 10134
    VVVV1_0( w_fp[662], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 10135 OF 15495 ***
    // Wavefunction(s) for diagram number 10135
    // (none)
    // Amplitude(s) for diagram number 10135
    VVV1_0( w_fp[150], w_fp[7], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10136 OF 15495 ***
    // Wavefunction(s) for diagram number 10136
    // (none)
    // Amplitude(s) for diagram number 10136
    VVV1_0( w_fp[150], w_fp[5], w_fp[680], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 10137 OF 15495 ***
    // Wavefunction(s) for diagram number 10137
    // (none)
    // Amplitude(s) for diagram number 10137
    FFV1_0( w_fp[175], w_fp[493], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];

    // *** DIAGRAM 10138 OF 15495 ***
    // Wavefunction(s) for diagram number 10138
    // (none)
    // Amplitude(s) for diagram number 10138
    FFV1_0( w_fp[175], w_fp[2], w_fp[680], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10139 OF 15495 ***
    // Wavefunction(s) for diagram number 10139
    // (none)
    // Amplitude(s) for diagram number 10139
    FFV1_0( w_fp[176], w_fp[493], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];

    // *** DIAGRAM 10140 OF 15495 ***
    // Wavefunction(s) for diagram number 10140
    // (none)
    // Amplitude(s) for diagram number 10140
    FFV1_0( w_fp[176], w_fp[2], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10141 OF 15495 ***
    // Wavefunction(s) for diagram number 10141
    // (none)
    // Amplitude(s) for diagram number 10141
    FFV1_0( w_fp[577], w_fp[483], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10142 OF 15495 ***
    // Wavefunction(s) for diagram number 10142
    // (none)
    // Amplitude(s) for diagram number 10142
    FFV1_0( w_fp[578], w_fp[483], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10143 OF 15495 ***
    // Wavefunction(s) for diagram number 10143
    FFV1P0_3( w_fp[594], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[724] );
    // Amplitude(s) for diagram number 10143
    VVV1_0( w_fp[260], w_fp[7], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10144 OF 15495 ***
    // Wavefunction(s) for diagram number 10144
    // (none)
    // Amplitude(s) for diagram number 10144
    FFV1_0( w_fp[578], w_fp[2], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];

    // *** DIAGRAM 10145 OF 15495 ***
    // Wavefunction(s) for diagram number 10145
    // (none)
    // Amplitude(s) for diagram number 10145
    VVV1_0( w_fp[265], w_fp[5], w_fp[724], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10146 OF 15495 ***
    // Wavefunction(s) for diagram number 10146
    // (none)
    // Amplitude(s) for diagram number 10146
    FFV1_0( w_fp[577], w_fp[2], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];

    // *** DIAGRAM 10147 OF 15495 ***
    // Wavefunction(s) for diagram number 10147
    // (none)
    // Amplitude(s) for diagram number 10147
    FFV1_0( w_fp[594], w_fp[2], w_fp[269], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[270], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[594], w_fp[2], w_fp[271], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10148 OF 15495 ***
    // Wavefunction(s) for diagram number 10148
    // (none)
    // Amplitude(s) for diagram number 10148
    FFV1_0( w_fp[175], w_fp[748], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10149 OF 15495 ***
    // Wavefunction(s) for diagram number 10149
    // (none)
    // Amplitude(s) for diagram number 10149
    FFV1_0( w_fp[530], w_fp[483], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10150 OF 15495 ***
    // Wavefunction(s) for diagram number 10150
    // (none)
    // Amplitude(s) for diagram number 10150
    FFV1_0( w_fp[176], w_fp[748], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10151 OF 15495 ***
    // Wavefunction(s) for diagram number 10151
    // (none)
    // Amplitude(s) for diagram number 10151
    FFV1_0( w_fp[555], w_fp[483], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10152 OF 15495 ***
    // Wavefunction(s) for diagram number 10152
    // (none)
    // Amplitude(s) for diagram number 10152
    VVVV1_0( w_fp[0], w_fp[260], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[260], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[260], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];

    // *** DIAGRAM 10153 OF 15495 ***
    // Wavefunction(s) for diagram number 10153
    // (none)
    // Amplitude(s) for diagram number 10153
    VVV1_0( w_fp[150], w_fp[7], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10154 OF 15495 ***
    // Wavefunction(s) for diagram number 10154
    VVV1P0_1( w_fp[0], w_fp[150], COUPs[0], 1.0, depCoup, 0., 0., w_fp[725] );
    // Amplitude(s) for diagram number 10154
    VVV1_0( w_fp[260], w_fp[7], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10155 OF 15495 ***
    // Wavefunction(s) for diagram number 10155
    // (none)
    // Amplitude(s) for diagram number 10155
    FFV1_0( w_fp[176], w_fp[2], w_fp[700], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10156 OF 15495 ***
    // Wavefunction(s) for diagram number 10156
    // (none)
    // Amplitude(s) for diagram number 10156
    FFV1_0( w_fp[555], w_fp[2], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];

    // *** DIAGRAM 10157 OF 15495 ***
    // Wavefunction(s) for diagram number 10157
    // (none)
    // Amplitude(s) for diagram number 10157
    VVVV1_0( w_fp[0], w_fp[265], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[265], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[265], w_fp[150], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 10158 OF 15495 ***
    // Wavefunction(s) for diagram number 10158
    // (none)
    // Amplitude(s) for diagram number 10158
    VVV1_0( w_fp[150], w_fp[5], w_fp[705], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];

    // *** DIAGRAM 10159 OF 15495 ***
    // Wavefunction(s) for diagram number 10159
    // (none)
    // Amplitude(s) for diagram number 10159
    VVV1_0( w_fp[265], w_fp[5], w_fp[725], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];

    // *** DIAGRAM 10160 OF 15495 ***
    // Wavefunction(s) for diagram number 10160
    // (none)
    // Amplitude(s) for diagram number 10160
    FFV1_0( w_fp[175], w_fp[2], w_fp[705], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10161 OF 15495 ***
    // Wavefunction(s) for diagram number 10161
    // (none)
    // Amplitude(s) for diagram number 10161
    FFV1_0( w_fp[530], w_fp[2], w_fp[265], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 10162 OF 15495 ***
    // Wavefunction(s) for diagram number 10162
    // (none)
    // Amplitude(s) for diagram number 10162
    VVV1_0( w_fp[706], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVV1_0( w_fp[707], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    VVV1_0( w_fp[708], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 10163 OF 15495 ***
    // Wavefunction(s) for diagram number 10163
    // (none)
    // Amplitude(s) for diagram number 10163
    FFV1_0( w_fp[176], w_fp[2], w_fp[706], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[707], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[708], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10164 OF 15495 ***
    // Wavefunction(s) for diagram number 10164
    // (none)
    // Amplitude(s) for diagram number 10164
    VVV1_0( w_fp[712], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    VVV1_0( w_fp[713], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    VVV1_0( w_fp[714], w_fp[150], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[600] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];

    // *** DIAGRAM 10165 OF 15495 ***
    // Wavefunction(s) for diagram number 10165
    // (none)
    // Amplitude(s) for diagram number 10165
    FFV1_0( w_fp[175], w_fp[2], w_fp[712], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[713], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[175], w_fp[2], w_fp[714], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[600] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10166 OF 15495 ***
    // Wavefunction(s) for diagram number 10166
    // (none)
    // Amplitude(s) for diagram number 10166
    VVV1_0( w_fp[0], w_fp[269], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[270], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[271], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10167 OF 15495 ***
    // Wavefunction(s) for diagram number 10167
    // (none)
    // Amplitude(s) for diagram number 10167
    VVV1_0( w_fp[662], w_fp[150], w_fp[100], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10168 OF 15495 ***
    // Wavefunction(s) for diagram number 10168
    // (none)
    // Amplitude(s) for diagram number 10168
    FFV1_0( w_fp[174], w_fp[122], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10169 OF 15495 ***
    // Wavefunction(s) for diagram number 10169
    // (none)
    // Amplitude(s) for diagram number 10169
    FFV1_0( w_fp[109], w_fp[2], w_fp[662], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[4] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10170 OF 15495 ***
    // Wavefunction(s) for diagram number 10170
    // (none)
    // Amplitude(s) for diagram number 10170
    FFV1_0( w_fp[594], w_fp[483], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];

    // *** DIAGRAM 10171 OF 15495 ***
    // Wavefunction(s) for diagram number 10171
    // (none)
    // Amplitude(s) for diagram number 10171
    FFV1_0( w_fp[594], w_fp[2], w_fp[287], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10172 OF 15495 ***
    // Wavefunction(s) for diagram number 10172
    // (none)
    // Amplitude(s) for diagram number 10172
    FFV1_0( w_fp[594], w_fp[122], w_fp[258], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10173 OF 15495 ***
    // Wavefunction(s) for diagram number 10173
    // (none)
    // Amplitude(s) for diagram number 10173
    FFV1_0( w_fp[174], w_fp[483], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[268] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10174 OF 15495 ***
    // Wavefunction(s) for diagram number 10174
    // (none)
    // Amplitude(s) for diagram number 10174
    FFV1_0( w_fp[283], w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10175 OF 15495 ***
    // Wavefunction(s) for diagram number 10175
    // (none)
    // Amplitude(s) for diagram number 10175
    VVV1_0( w_fp[522], w_fp[258], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 10176 OF 15495 ***
    // Wavefunction(s) for diagram number 10176
    // (none)
    // Amplitude(s) for diagram number 10176
    FFV1_0( w_fp[109], w_fp[483], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];

    // *** DIAGRAM 10177 OF 15495 ***
    // Wavefunction(s) for diagram number 10177
    // (none)
    // Amplitude(s) for diagram number 10177
    FFV1_0( w_fp[283], w_fp[122], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 10178 OF 15495 ***
    // Wavefunction(s) for diagram number 10178
    // (none)
    // Amplitude(s) for diagram number 10178
    VVV1_0( w_fp[0], w_fp[287], w_fp[150], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 10179 OF 15495 ***
    // Wavefunction(s) for diagram number 10179
    // (none)
    // Amplitude(s) for diagram number 10179
    FFV1_0( w_fp[174], w_fp[2], w_fp[733], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[734], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    FFV1_0( w_fp[174], w_fp[2], w_fp[735], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 10180 OF 15495 ***
    // Wavefunction(s) for diagram number 10180
    // (none)
    // Amplitude(s) for diagram number 10180
    VVVV1_0( w_fp[662], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[2] -= amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[26] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[266] -= amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV3_0( w_fp[662], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0( w_fp[662], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 10181 OF 15495 ***
    // Wavefunction(s) for diagram number 10181
    // (none)
    // Amplitude(s) for diagram number 10181
    VVV1_0( w_fp[144], w_fp[6], w_fp[544], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    jamp_sv[552] += amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 10182 OF 15495 ***
    // Wavefunction(s) for diagram number 10182
    // (none)
    // Amplitude(s) for diagram number 10182
    VVV1_0( w_fp[144], w_fp[5], w_fp[670], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[366] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[434] -= amp_sv[0];
    jamp_sv[440] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[510] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 10183 OF 15495 ***
    // Wavefunction(s) for diagram number 10183
    // (none)
    // Amplitude(s) for diagram number 10183
    FFV1_0( w_fp[180], w_fp[493], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[26] -= amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[266] += amp_sv[0];

    // *** DIAGRAM 10184 OF 15495 ***
    // Wavefunction(s) for diagram number 10184
    // (none)
    // Amplitude(s) for diagram number 10184
    FFV1_0( w_fp[180], w_fp[2], w_fp[670], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10185 OF 15495 ***
    // Wavefunction(s) for diagram number 10185
    // (none)
    // Amplitude(s) for diagram number 10185
    FFV1_0( w_fp[181], w_fp[493], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[264] += amp_sv[0];

    // *** DIAGRAM 10186 OF 15495 ***
    // Wavefunction(s) for diagram number 10186
    // (none)
    // Amplitude(s) for diagram number 10186
    FFV1_0( w_fp[181], w_fp[2], w_fp[544], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10187 OF 15495 ***
    // Wavefunction(s) for diagram number 10187
    // (none)
    // Amplitude(s) for diagram number 10187
    FFV1_0( w_fp[438], w_fp[483], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10188 OF 15495 ***
    // Wavefunction(s) for diagram number 10188
    // (none)
    // Amplitude(s) for diagram number 10188
    FFV1_0( w_fp[495], w_fp[483], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10189 OF 15495 ***
    // Wavefunction(s) for diagram number 10189
    FFV1P0_3( w_fp[596], w_fp[2], COUPs[1], 1.0, depCoup, 0., 0., w_fp[735] );
    // Amplitude(s) for diagram number 10189
    VVV1_0( w_fp[260], w_fp[6], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10190 OF 15495 ***
    // Wavefunction(s) for diagram number 10190
    // (none)
    // Amplitude(s) for diagram number 10190
    FFV1_0( w_fp[495], w_fp[2], w_fp[260], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];

    // *** DIAGRAM 10191 OF 15495 ***
    // Wavefunction(s) for diagram number 10191
    // (none)
    // Amplitude(s) for diagram number 10191
    VVV1_0( w_fp[263], w_fp[5], w_fp[735], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10192 OF 15495 ***
    // Wavefunction(s) for diagram number 10192
    // (none)
    // Amplitude(s) for diagram number 10192
    FFV1_0( w_fp[438], w_fp[2], w_fp[263], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += amp_sv[0];
    jamp_sv[278] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];

    // *** DIAGRAM 10193 OF 15495 ***
    // Wavefunction(s) for diagram number 10193
    // (none)
    // Amplitude(s) for diagram number 10193
    FFV1_0( w_fp[596], w_fp[2], w_fp[266], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[267], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[278] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[440] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[596], w_fp[2], w_fp[268], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10194 OF 15495 ***
    // Wavefunction(s) for diagram number 10194
    // (none)
    // Amplitude(s) for diagram number 10194
    FFV1_0( w_fp[180], w_fp[748], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[266] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10195 OF 15495 ***
    // Wavefunction(s) for diagram number 10195
    // (none)
    // Amplitude(s) for diagram number 10195
    FFV1_0( w_fp[557], w_fp[483], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10196 OF 15495 ***
    // Wavefunction(s) for diagram number 10196
    // (none)
    // Amplitude(s) for diagram number 10196
    FFV1_0( w_fp[181], w_fp[748], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[264] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10197 OF 15495 ***
    // Wavefunction(s) for diagram number 10197
    // (none)
    // Amplitude(s) for diagram number 10197
    FFV1_0( w_fp[583], w_fp[483], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10198 OF 15495 ***
    // Wavefunction(s) for diagram number 10198
    // (none)
    // Amplitude(s) for diagram number 10198
    VVVV1_0( w_fp[0], w_fp[260], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[260], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[260], w_fp[144], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];

    // *** DIAGRAM 10199 OF 15495 ***
    // Wavefunction(s) for diagram number 10199
    // (none)
    // Amplitude(s) for diagram number 10199
    VVV1_0( w_fp[144], w_fp[6], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 10200 OF 15495 ***
    // Wavefunction(s) for diagram number 10200
    VVV1P0_1( w_fp[0], w_fp[144], COUPs[0], 1.0, depCoup, 0., 0., w_fp[734] );
    // Amplitude(s) for diagram number 10200
    VVV1_0( w_fp[260], w_fp[6], w_fp[734], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
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
    storeWf( wfs, w_cx, nevt, 491 );
    storeWf( wfs, w_cx, nevt, 493 );
    storeWf( wfs, w_cx, nevt, 724 );
    storeWf( wfs, w_cx, nevt, 725 );
    storeWf( wfs, w_cx, nevt, 734 );
    storeWf( wfs, w_cx, nevt, 735 );
    storeWf( wfs, w_cx, nevt, 746 );
    storeWf( wfs, w_cx, nevt, 747 );
    storeWf( wfs, w_cx, nevt, 748 );
    storeWf( wfs, w_cx, nevt, 749 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
