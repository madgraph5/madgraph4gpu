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
  diagramgroup85( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 170 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 214 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 243 );
    retrieveWf( wfs, w_cx, nevt, 244 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 359 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 526 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 593 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 596 );
#endif
#endif

    // *** DIAGRAM 8401 OF 15495 ***
    // Wavefunction(s) for diagram number 8401
    // (none)
    // Amplitude(s) for diagram number 8401
    VVV1_0( w_fp[594], w_fp[1], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[513] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[544] -= amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    jamp_sv[568] += amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[581] -= amp_sv[0];

    // *** DIAGRAM 8402 OF 15495 ***
    // Wavefunction(s) for diagram number 8402
    // (none)
    // Amplitude(s) for diagram number 8402
    FFV1_0( w_fp[221], w_fp[476], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];

    // *** DIAGRAM 8403 OF 15495 ***
    // Wavefunction(s) for diagram number 8403
    // (none)
    // Amplitude(s) for diagram number 8403
    VVV1_0( w_fp[532], w_fp[254], w_fp[214], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];

    // *** DIAGRAM 8404 OF 15495 ***
    // Wavefunction(s) for diagram number 8404
    // (none)
    // Amplitude(s) for diagram number 8404
    FFV1_0( w_fp[3], w_fp[197], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    jamp_sv[567] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[501] += amp_sv[0];
    jamp_sv[503] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[513] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[519] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[581] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[197], w_fp[593], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[544] += amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[568] -= amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];

    // *** DIAGRAM 8405 OF 15495 ***
    // Wavefunction(s) for diagram number 8405
    // (none)
    // Amplitude(s) for diagram number 8405
    VVVV1_0( w_fp[585], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];
    VVVV3_0( w_fp[585], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    VVVV4_0( w_fp[585], w_fp[239], w_fp[5], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];

    // *** DIAGRAM 8406 OF 15495 ***
    // Wavefunction(s) for diagram number 8406
    // (none)
    // Amplitude(s) for diagram number 8406
    VVV1_0( w_fp[239], w_fp[6], w_fp[576], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];

    // *** DIAGRAM 8407 OF 15495 ***
    // Wavefunction(s) for diagram number 8407
    // (none)
    // Amplitude(s) for diagram number 8407
    VVV1_0( w_fp[239], w_fp[5], w_fp[556], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];

    // *** DIAGRAM 8408 OF 15495 ***
    // Wavefunction(s) for diagram number 8408
    FFV1_1( w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[437] );
    // Amplitude(s) for diagram number 8408
    FFV1_0( w_fp[216], w_fp[437], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];

    // *** DIAGRAM 8409 OF 15495 ***
    // Wavefunction(s) for diagram number 8409
    // (none)
    // Amplitude(s) for diagram number 8409
    FFV1_0( w_fp[216], w_fp[2], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8410 OF 15495 ***
    // Wavefunction(s) for diagram number 8410
    // (none)
    // Amplitude(s) for diagram number 8410
    FFV1_0( w_fp[198], w_fp[437], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];

    // *** DIAGRAM 8411 OF 15495 ***
    // Wavefunction(s) for diagram number 8411
    // (none)
    // Amplitude(s) for diagram number 8411
    FFV1_0( w_fp[198], w_fp[2], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8412 OF 15495 ***
    // Wavefunction(s) for diagram number 8412
    // (none)
    // Amplitude(s) for diagram number 8412
    FFV1_0( w_fp[355], w_fp[481], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8413 OF 15495 ***
    // Wavefunction(s) for diagram number 8413
    // (none)
    // Amplitude(s) for diagram number 8413
    FFV1_0( w_fp[355], w_fp[579], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8414 OF 15495 ***
    // Wavefunction(s) for diagram number 8414
    FFV1_1( w_fp[453], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[571] );
    // Amplitude(s) for diagram number 8414
    FFV1_0( w_fp[216], w_fp[571], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8415 OF 15495 ***
    // Wavefunction(s) for diagram number 8415
    // (none)
    // Amplitude(s) for diagram number 8415
    FFV1_0( w_fp[216], w_fp[579], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8416 OF 15495 ***
    // Wavefunction(s) for diagram number 8416
    // (none)
    // Amplitude(s) for diagram number 8416
    FFV1_0( w_fp[198], w_fp[571], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8417 OF 15495 ***
    // Wavefunction(s) for diagram number 8417
    // (none)
    // Amplitude(s) for diagram number 8417
    FFV1_0( w_fp[198], w_fp[481], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8418 OF 15495 ***
    // Wavefunction(s) for diagram number 8418
    // (none)
    // Amplitude(s) for diagram number 8418
    FFV1_0( w_fp[355], w_fp[509], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];

    // *** DIAGRAM 8419 OF 15495 ***
    // Wavefunction(s) for diagram number 8419
    // (none)
    // Amplitude(s) for diagram number 8419
    FFV1_0( w_fp[355], w_fp[2], w_fp[557], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8420 OF 15495 ***
    // Wavefunction(s) for diagram number 8420
    // (none)
    // Amplitude(s) for diagram number 8420
    VVVV1_0( w_fp[505], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    VVVV3_0( w_fp[505], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += amp_sv[0];
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    VVVV4_0( w_fp[505], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];

    // *** DIAGRAM 8421 OF 15495 ***
    // Wavefunction(s) for diagram number 8421
    // (none)
    // Amplitude(s) for diagram number 8421
    VVV1_0( w_fp[239], w_fp[6], w_fp[520], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[520] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];

    // *** DIAGRAM 8422 OF 15495 ***
    // Wavefunction(s) for diagram number 8422
    // (none)
    // Amplitude(s) for diagram number 8422
    VVV1_0( w_fp[1], w_fp[239], w_fp[557], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];

    // *** DIAGRAM 8423 OF 15495 ***
    // Wavefunction(s) for diagram number 8423
    // (none)
    // Amplitude(s) for diagram number 8423
    FFV1_0( w_fp[198], w_fp[2], w_fp[520], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8424 OF 15495 ***
    // Wavefunction(s) for diagram number 8424
    // (none)
    // Amplitude(s) for diagram number 8424
    FFV1_0( w_fp[198], w_fp[509], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];

    // *** DIAGRAM 8425 OF 15495 ***
    // Wavefunction(s) for diagram number 8425
    // (none)
    // Amplitude(s) for diagram number 8425
    FFV1_0( w_fp[355], w_fp[518], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];

    // *** DIAGRAM 8426 OF 15495 ***
    // Wavefunction(s) for diagram number 8426
    // (none)
    // Amplitude(s) for diagram number 8426
    FFV1_0( w_fp[355], w_fp[2], w_fp[435], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8427 OF 15495 ***
    // Wavefunction(s) for diagram number 8427
    // (none)
    // Amplitude(s) for diagram number 8427
    VVVV1_0( w_fp[474], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[619] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[1], w_fp[239], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 8428 OF 15495 ***
    // Wavefunction(s) for diagram number 8428
    // (none)
    // Amplitude(s) for diagram number 8428
    VVV1_0( w_fp[239], w_fp[5], w_fp[486], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];

    // *** DIAGRAM 8429 OF 15495 ***
    // Wavefunction(s) for diagram number 8429
    // (none)
    // Amplitude(s) for diagram number 8429
    VVV1_0( w_fp[1], w_fp[239], w_fp[435], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 8430 OF 15495 ***
    // Wavefunction(s) for diagram number 8430
    // (none)
    // Amplitude(s) for diagram number 8430
    FFV1_0( w_fp[216], w_fp[2], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8431 OF 15495 ***
    // Wavefunction(s) for diagram number 8431
    // (none)
    // Amplitude(s) for diagram number 8431
    FFV1_0( w_fp[216], w_fp[518], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[619] -= amp_sv[0];

    // *** DIAGRAM 8432 OF 15495 ***
    // Wavefunction(s) for diagram number 8432
    // (none)
    // Amplitude(s) for diagram number 8432
    VVV1_0( w_fp[543], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    VVV1_0( w_fp[538], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    VVV1_0( w_fp[527], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];

    // *** DIAGRAM 8433 OF 15495 ***
    // Wavefunction(s) for diagram number 8433
    // (none)
    // Amplitude(s) for diagram number 8433
    FFV1_0( w_fp[198], w_fp[2], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[109] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8434 OF 15495 ***
    // Wavefunction(s) for diagram number 8434
    // (none)
    // Amplitude(s) for diagram number 8434
    VVV1_0( w_fp[558], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[509] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    VVV1_0( w_fp[480], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[460] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[619] += amp_sv[0];
    VVV1_0( w_fp[550], w_fp[239], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 8435 OF 15495 ***
    // Wavefunction(s) for diagram number 8435
    // (none)
    // Amplitude(s) for diagram number 8435
    FFV1_0( w_fp[216], w_fp[2], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[216], w_fp[2], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8436 OF 15495 ***
    // Wavefunction(s) for diagram number 8436
    // (none)
    // Amplitude(s) for diagram number 8436
    FFV1_0( w_fp[355], w_fp[2], w_fp[495], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[438], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[596], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8437 OF 15495 ***
    // Wavefunction(s) for diagram number 8437
    // (none)
    // Amplitude(s) for diagram number 8437
    VVV1_0( w_fp[495], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    VVV1_0( w_fp[438], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[118] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];
    VVV1_0( w_fp[596], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[460] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];

    // *** DIAGRAM 8438 OF 15495 ***
    // Wavefunction(s) for diagram number 8438
    // (none)
    // Amplitude(s) for diagram number 8438
    VVV1_0( w_fp[585], w_fp[239], w_fp[113], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 8439 OF 15495 ***
    // Wavefunction(s) for diagram number 8439
    // (none)
    // Amplitude(s) for diagram number 8439
    FFV1_0( w_fp[196], w_fp[244], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8440 OF 15495 ***
    // Wavefunction(s) for diagram number 8440
    // (none)
    // Amplitude(s) for diagram number 8440
    FFV1_0( w_fp[243], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8441 OF 15495 ***
    // Wavefunction(s) for diagram number 8441
    // (none)
    // Amplitude(s) for diagram number 8441
    FFV1_0( w_fp[355], w_fp[453], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 8442 OF 15495 ***
    // Wavefunction(s) for diagram number 8442
    // (none)
    // Amplitude(s) for diagram number 8442
    FFV1_0( w_fp[196], w_fp[453], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[112] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[603] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8443 OF 15495 ***
    // Wavefunction(s) for diagram number 8443
    // (none)
    // Amplitude(s) for diagram number 8443
    FFV1_0( w_fp[243], w_fp[453], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 8444 OF 15495 ***
    // Wavefunction(s) for diagram number 8444
    // (none)
    // Amplitude(s) for diagram number 8444
    FFV1_0( w_fp[442], w_fp[2], w_fp[359], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[181] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[208] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[442] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[562] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8445 OF 15495 ***
    // Wavefunction(s) for diagram number 8445
    // (none)
    // Amplitude(s) for diagram number 8445
    FFV1_0( w_fp[442], w_fp[244], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];

    // *** DIAGRAM 8446 OF 15495 ***
    // Wavefunction(s) for diagram number 8446
    // (none)
    // Amplitude(s) for diagram number 8446
    FFV1_0( w_fp[355], w_fp[2], w_fp[531], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[436] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[622] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8447 OF 15495 ***
    // Wavefunction(s) for diagram number 8447
    // (none)
    // Amplitude(s) for diagram number 8447
    VVV1_0( w_fp[531], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[208] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[622] -= amp_sv[0];

    // *** DIAGRAM 8448 OF 15495 ***
    // Wavefunction(s) for diagram number 8448
    // (none)
    // Amplitude(s) for diagram number 8448
    FFV1_0( w_fp[355], w_fp[244], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[436] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[570] += amp_sv[0];

    // *** DIAGRAM 8449 OF 15495 ***
    // Wavefunction(s) for diagram number 8449
    // (none)
    // Amplitude(s) for diagram number 8449
    VVV1_0( w_fp[532], w_fp[359], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];

    // *** DIAGRAM 8450 OF 15495 ***
    // Wavefunction(s) for diagram number 8450
    // (none)
    // Amplitude(s) for diagram number 8450
    FFV1_0( w_fp[196], w_fp[2], w_fp[539], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[442] += amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[562] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[605] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[591], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[112] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[208] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[622] += amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[592], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[99] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[141] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[436] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[442] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    jamp_sv[562] += amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[605] += amp_sv[0];

    // *** DIAGRAM 8451 OF 15495 ***
    // Wavefunction(s) for diagram number 8451
    // (none)
    // Amplitude(s) for diagram number 8451
    VVVV1_0( w_fp[585], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];
    VVVV3_0( w_fp[585], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    VVVV4_0( w_fp[585], w_fp[154], w_fp[4], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];

    // *** DIAGRAM 8452 OF 15495 ***
    // Wavefunction(s) for diagram number 8452
    // (none)
    // Amplitude(s) for diagram number 8452
    VVV1_0( w_fp[154], w_fp[6], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];

    // *** DIAGRAM 8453 OF 15495 ***
    // Wavefunction(s) for diagram number 8453
    // (none)
    // Amplitude(s) for diagram number 8453
    VVV1_0( w_fp[154], w_fp[4], w_fp[556], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];

    // *** DIAGRAM 8454 OF 15495 ***
    // Wavefunction(s) for diagram number 8454
    // (none)
    // Amplitude(s) for diagram number 8454
    FFV1_0( w_fp[218], w_fp[437], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];

    // *** DIAGRAM 8455 OF 15495 ***
    // Wavefunction(s) for diagram number 8455
    // (none)
    // Amplitude(s) for diagram number 8455
    FFV1_0( w_fp[218], w_fp[2], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8456 OF 15495 ***
    // Wavefunction(s) for diagram number 8456
    // (none)
    // Amplitude(s) for diagram number 8456
    FFV1_0( w_fp[170], w_fp[437], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[601] -= amp_sv[0];

    // *** DIAGRAM 8457 OF 15495 ***
    // Wavefunction(s) for diagram number 8457
    // (none)
    // Amplitude(s) for diagram number 8457
    FFV1_0( w_fp[170], w_fp[2], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8458 OF 15495 ***
    // Wavefunction(s) for diagram number 8458
    // (none)
    // Amplitude(s) for diagram number 8458
    FFV1_0( w_fp[256], w_fp[452], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8459 OF 15495 ***
    // Wavefunction(s) for diagram number 8459
    // (none)
    // Amplitude(s) for diagram number 8459
    FFV1_0( w_fp[256], w_fp[579], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8460 OF 15495 ***
    // Wavefunction(s) for diagram number 8460
    // (none)
    // Amplitude(s) for diagram number 8460
    FFV1_0( w_fp[218], w_fp[571], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8461 OF 15495 ***
    // Wavefunction(s) for diagram number 8461
    // (none)
    // Amplitude(s) for diagram number 8461
    FFV1_0( w_fp[218], w_fp[579], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8462 OF 15495 ***
    // Wavefunction(s) for diagram number 8462
    // (none)
    // Amplitude(s) for diagram number 8462
    FFV1_0( w_fp[170], w_fp[571], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8463 OF 15495 ***
    // Wavefunction(s) for diagram number 8463
    // (none)
    // Amplitude(s) for diagram number 8463
    FFV1_0( w_fp[170], w_fp[452], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8464 OF 15495 ***
    // Wavefunction(s) for diagram number 8464
    // (none)
    // Amplitude(s) for diagram number 8464
    FFV1_0( w_fp[256], w_fp[444], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];

    // *** DIAGRAM 8465 OF 15495 ***
    // Wavefunction(s) for diagram number 8465
    // (none)
    // Amplitude(s) for diagram number 8465
    FFV1_0( w_fp[256], w_fp[2], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8466 OF 15495 ***
    // Wavefunction(s) for diagram number 8466
    // (none)
    // Amplitude(s) for diagram number 8466
    VVVV1_0( w_fp[537], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    VVVV3_0( w_fp[537], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    VVVV4_0( w_fp[537], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];

    // *** DIAGRAM 8467 OF 15495 ***
    // Wavefunction(s) for diagram number 8467
    // (none)
    // Amplitude(s) for diagram number 8467
    VVV1_0( w_fp[154], w_fp[6], w_fp[522], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[514] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];

    // *** DIAGRAM 8468 OF 15495 ***
    // Wavefunction(s) for diagram number 8468
    // (none)
    // Amplitude(s) for diagram number 8468
    VVV1_0( w_fp[1], w_fp[154], w_fp[553], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];

    // *** DIAGRAM 8469 OF 15495 ***
    // Wavefunction(s) for diagram number 8469
    // (none)
    // Amplitude(s) for diagram number 8469
    FFV1_0( w_fp[170], w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8470 OF 15495 ***
    // Wavefunction(s) for diagram number 8470
    // (none)
    // Amplitude(s) for diagram number 8470
    FFV1_0( w_fp[170], w_fp[444], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];

    // *** DIAGRAM 8471 OF 15495 ***
    // Wavefunction(s) for diagram number 8471
    // (none)
    // Amplitude(s) for diagram number 8471
    FFV1_0( w_fp[256], w_fp[518], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];

    // *** DIAGRAM 8472 OF 15495 ***
    // Wavefunction(s) for diagram number 8472
    // (none)
    // Amplitude(s) for diagram number 8472
    FFV1_0( w_fp[256], w_fp[2], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8473 OF 15495 ***
    // Wavefunction(s) for diagram number 8473
    // (none)
    // Amplitude(s) for diagram number 8473
    VVVV1_0( w_fp[474], w_fp[1], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    VVVV3_0( w_fp[474], w_fp[1], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[618] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    VVVV4_0( w_fp[474], w_fp[1], w_fp[154], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 8474 OF 15495 ***
    // Wavefunction(s) for diagram number 8474
    // (none)
    // Amplitude(s) for diagram number 8474
    VVV1_0( w_fp[154], w_fp[4], w_fp[486], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[277] -= amp_sv[0];
    jamp_sv[280] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];

    // *** DIAGRAM 8475 OF 15495 ***
    // Wavefunction(s) for diagram number 8475
    // (none)
    // Amplitude(s) for diagram number 8475
    VVV1_0( w_fp[1], w_fp[154], w_fp[582], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 8476 OF 15495 ***
    // Wavefunction(s) for diagram number 8476
    // (none)
    // Amplitude(s) for diagram number 8476
    FFV1_0( w_fp[218], w_fp[2], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8477 OF 15495 ***
    // Wavefunction(s) for diagram number 8477
    // (none)
    // Amplitude(s) for diagram number 8477
    FFV1_0( w_fp[218], w_fp[518], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[618] -= amp_sv[0];

    // *** DIAGRAM 8478 OF 15495 ***
    // Wavefunction(s) for diagram number 8478
    // (none)
    // Amplitude(s) for diagram number 8478
    VVV1_0( w_fp[536], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[269] += amp_sv[0];
    jamp_sv[283] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[576] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    VVV1_0( w_fp[560], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[163] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    VVV1_0( w_fp[559], w_fp[154], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];

    // *** DIAGRAM 8479 OF 15495 ***
    // Wavefunction(s) for diagram number 8479
    // (none)
    // Amplitude(s) for diagram number 8479
    FFV1_0( w_fp[170], w_fp[2], w_fp[536], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[103] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[163] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[170], w_fp[2], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[269] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[283] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8480 OF 15495 ***
    // Wavefunction(s) for diagram number 8480
    // (none)
    // Amplitude(s) for diagram number 8480
    VVV1_0( w_fp[558], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += amp_sv[0];
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[337] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[508] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    VVV1_0( w_fp[480], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[277] += amp_sv[0];
    jamp_sv[280] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[618] += amp_sv[0];
    VVV1_0( w_fp[550], w_fp[154], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[269] -= amp_sv[0];
    jamp_sv[283] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[337] -= amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[576] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];

    // *** DIAGRAM 8481 OF 15495 ***
    // Wavefunction(s) for diagram number 8481
    // (none)
    // Amplitude(s) for diagram number 8481
    FFV1_0( w_fp[218], w_fp[2], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[114] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[218], w_fp[2], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[576] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8482 OF 15495 ***
    // Wavefunction(s) for diagram number 8482
    // (none)
    // Amplitude(s) for diagram number 8482
    FFV1_0( w_fp[256], w_fp[2], w_fp[471], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[256], w_fp[2], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8483 OF 15495 ***
    // Wavefunction(s) for diagram number 8483
    // (none)
    // Amplitude(s) for diagram number 8483
    VVV1_0( w_fp[471], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    VVV1_0( w_fp[562], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[116] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    VVV1_0( w_fp[524], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[163] -= amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];

    // *** DIAGRAM 8484 OF 15495 ***
    // Wavefunction(s) for diagram number 8484
    // (none)
    // Amplitude(s) for diagram number 8484
    VVV1_0( w_fp[585], w_fp[154], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];

    // *** DIAGRAM 8485 OF 15495 ***
    // Wavefunction(s) for diagram number 8485
    // (none)
    // Amplitude(s) for diagram number 8485
    FFV1_0( w_fp[168], w_fp[118], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8486 OF 15495 ***
    // Wavefunction(s) for diagram number 8486
    // (none)
    // Amplitude(s) for diagram number 8486
    FFV1_0( w_fp[200], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8487 OF 15495 ***
    // Wavefunction(s) for diagram number 8487
    // (none)
    // Amplitude(s) for diagram number 8487
    FFV1_0( w_fp[256], w_fp[453], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 8488 OF 15495 ***
    // Wavefunction(s) for diagram number 8488
    // (none)
    // Amplitude(s) for diagram number 8488
    FFV1_0( w_fp[168], w_fp[453], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[601] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8489 OF 15495 ***
    // Wavefunction(s) for diagram number 8489
    // (none)
    // Amplitude(s) for diagram number 8489
    FFV1_0( w_fp[200], w_fp[453], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];

    // *** DIAGRAM 8490 OF 15495 ***
    // Wavefunction(s) for diagram number 8490
    // (none)
    // Amplitude(s) for diagram number 8490
    FFV1_0( w_fp[45], w_fp[2], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[157] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[535] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8491 OF 15495 ***
    // Wavefunction(s) for diagram number 8491
    // (none)
    // Amplitude(s) for diagram number 8491
    FFV1_0( w_fp[45], w_fp[118], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];

    // *** DIAGRAM 8492 OF 15495 ***
    // Wavefunction(s) for diagram number 8492
    // (none)
    // Amplitude(s) for diagram number 8492
    FFV1_0( w_fp[256], w_fp[2], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[620] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8493 OF 15495 ***
    // Wavefunction(s) for diagram number 8493
    // (none)
    // Amplitude(s) for diagram number 8493
    VVV1_0( w_fp[526], w_fp[1], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[160] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[202] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[620] -= amp_sv[0];

    // *** DIAGRAM 8494 OF 15495 ***
    // Wavefunction(s) for diagram number 8494
    // (none)
    // Amplitude(s) for diagram number 8494
    FFV1_0( w_fp[256], w_fp[118], w_fp[532], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];

    // *** DIAGRAM 8495 OF 15495 ***
    // Wavefunction(s) for diagram number 8495
    // (none)
    // Amplitude(s) for diagram number 8495
    VVV1_0( w_fp[532], w_fp[362], w_fp[154], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];

    // *** DIAGRAM 8496 OF 15495 ***
    // Wavefunction(s) for diagram number 8496
    // (none)
    // Amplitude(s) for diagram number 8496
    FFV1_0( w_fp[168], w_fp[2], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[319] -= amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[535] += amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[601] += amp_sv[0];
    jamp_sv[604] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[2], w_fp[545], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[139] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[319] += amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[535] -= amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[601] -= amp_sv[0];
    jamp_sv[604] += amp_sv[0];

    // *** DIAGRAM 8497 OF 15495 ***
    // Wavefunction(s) for diagram number 8497
    // (none)
    // Amplitude(s) for diagram number 8497
    VVVV1_0( w_fp[585], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    VVVV3_0( w_fp[585], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];
    VVVV4_0( w_fp[585], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];

    // *** DIAGRAM 8498 OF 15495 ***
    // Wavefunction(s) for diagram number 8498
    // (none)
    // Amplitude(s) for diagram number 8498
    VVV1_0( w_fp[150], w_fp[5], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[96] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[600] -= amp_sv[0];

    // *** DIAGRAM 8499 OF 15495 ***
    // Wavefunction(s) for diagram number 8499
    // (none)
    // Amplitude(s) for diagram number 8499
    VVV1_0( w_fp[150], w_fp[4], w_fp[576], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];

    // *** DIAGRAM 8500 OF 15495 ***
    // Wavefunction(s) for diagram number 8500
    // (none)
    // Amplitude(s) for diagram number 8500
    FFV1_0( w_fp[202], w_fp[437], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 437 );
    storeWf( wfs, w_cx, nevt, 571 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
