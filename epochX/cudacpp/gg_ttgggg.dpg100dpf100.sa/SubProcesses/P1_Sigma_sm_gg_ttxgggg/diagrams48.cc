// Copyright (C) 2020-2025 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Sep 2025) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2025) for the MG5aMC CUDACPP plugin.

#include "GpuRuntime.h"
#include "HelAmps_sm.h"
#include "MemoryAccessAmplitudes.h"
#include "MemoryAccessChannelIds.h"
#include "MemoryAccessCouplings.h"
#include "MemoryAccessCouplingsFixed.h"
#include "MemoryAccessWavefunctions.h"
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
  diagramgroup48( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
#ifdef MGONGPUCPP_GPUIMPL
                  fptype* jamps,                  // output jamps[ncolor*2*nevt]
                  const int nGoodHel,             // input: number of good helicities
                  const fptype* couplings,        // input: dependent couplings[nevt*ndcoup*2] for all events
#else
                  cxtype_sv* jamps,               // output jamps[ncolor*2*neppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 112 );
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
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 437 );
    retrieveWf( wfs, w_cx, nevt, 444 );
    retrieveWf( wfs, w_cx, nevt, 445 );
    retrieveWf( wfs, w_cx, nevt, 448 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 451 );
    retrieveWf( wfs, w_cx, nevt, 452 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 471 );
    retrieveWf( wfs, w_cx, nevt, 472 );
    retrieveWf( wfs, w_cx, nevt, 474 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 490 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 535 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 549 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 561 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 582 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 587 );
#endif
#endif

    // *** DIAGRAM 4701 OF 15495 ***
    // Wavefunction(s) for diagram number 4701
    // (none)
    // Amplitude(s) for diagram number 4701
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[7], w_fp[571], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4702 OF 15495 ***
    // Wavefunction(s) for diagram number 4702
    // (none)
    // Amplitude(s) for diagram number 4702
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[154], w_fp[45], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4703 OF 15495 ***
    // Wavefunction(s) for diagram number 4703
    // (none)
    // Amplitude(s) for diagram number 4703
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[571], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4704 OF 15495 ***
    // Wavefunction(s) for diagram number 4704
    // (none)
    // Amplitude(s) for diagram number 4704
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[536], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[253] -= amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];

    // *** DIAGRAM 4705 OF 15495 ***
    // Wavefunction(s) for diagram number 4705
    // (none)
    // Amplitude(s) for diagram number 4705
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[537], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];

    // *** DIAGRAM 4706 OF 15495 ***
    // Wavefunction(s) for diagram number 4706
    // (none)
    // Amplitude(s) for diagram number 4706
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[509], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4707 OF 15495 ***
    // Wavefunction(s) for diagram number 4707
    // (none)
    // Amplitude(s) for diagram number 4707
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[154], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 4708 OF 15495 ***
    // Wavefunction(s) for diagram number 4708
    // (none)
    // Amplitude(s) for diagram number 4708
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[154], w_fp[6], w_fp[575], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4709 OF 15495 ***
    // Wavefunction(s) for diagram number 4709
    // (none)
    // Amplitude(s) for diagram number 4709
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[154], w_fp[509], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];

    // *** DIAGRAM 4710 OF 15495 ***
    // Wavefunction(s) for diagram number 4710
    // (none)
    // Amplitude(s) for diagram number 4710
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4711 OF 15495 ***
    // Wavefunction(s) for diagram number 4711
    // (none)
    // Amplitude(s) for diagram number 4711
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[537], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += amp_sv[0];
    jamp_sv[259] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];

    // *** DIAGRAM 4712 OF 15495 ***
    // Wavefunction(s) for diagram number 4712
    // (none)
    // Amplitude(s) for diagram number 4712
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[579], w_fp[154], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[505] += amp_sv[0];
    jamp_sv[511] -= amp_sv[0];
    jamp_sv[607] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[649] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[580], w_fp[154], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[193] -= amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[253] += amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[581], w_fp[154], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4713 OF 15495 ***
    // Wavefunction(s) for diagram number 4713
    // (none)
    // Amplitude(s) for diagram number 4713
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[579], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[580], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[199] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[2], w_fp[581], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[487] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[511] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4714 OF 15495 ***
    // Wavefunction(s) for diagram number 4714
    // (none)
    // Amplitude(s) for diagram number 4714
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[582], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += amp_sv[0];
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[487] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[631] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[223] += amp_sv[0];
    jamp_sv[259] += amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[154], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[487] += amp_sv[0];
    jamp_sv[505] -= amp_sv[0];
    jamp_sv[511] += amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[607] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[631] += amp_sv[0];
    jamp_sv[649] -= amp_sv[0];

    // *** DIAGRAM 4715 OF 15495 ***
    // Wavefunction(s) for diagram number 4715
    // (none)
    // Amplitude(s) for diagram number 4715
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[582], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[583], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[43] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[217] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[223] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[259] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[170], w_fp[2], w_fp[584], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[607] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[631] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4716 OF 15495 ***
    // Wavefunction(s) for diagram number 4716
    // (none)
    // Amplitude(s) for diagram number 4716
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[444], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[452], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[445], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[610] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4717 OF 15495 ***
    // Wavefunction(s) for diagram number 4717
    // (none)
    // Amplitude(s) for diagram number 4717
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[444], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[452], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[445], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[193] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[223] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4718 OF 15495 ***
    // Wavefunction(s) for diagram number 4718
    // (none)
    // Amplitude(s) for diagram number 4718
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[154], w_fp[84], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4719 OF 15495 ***
    // Wavefunction(s) for diagram number 4719
    // (none)
    // Amplitude(s) for diagram number 4719
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[128], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4720 OF 15495 ***
    // Wavefunction(s) for diagram number 4720
    // (none)
    // Amplitude(s) for diagram number 4720
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[112], w_fp[2], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4721 OF 15495 ***
    // Wavefunction(s) for diagram number 4721
    // (none)
    // Amplitude(s) for diagram number 4721
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[530], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];

    // *** DIAGRAM 4722 OF 15495 ***
    // Wavefunction(s) for diagram number 4722
    // (none)
    // Amplitude(s) for diagram number 4722
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[530], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[29] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[245] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4723 OF 15495 ***
    // Wavefunction(s) for diagram number 4723
    // (none)
    // Amplitude(s) for diagram number 4723
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[112], w_fp[530], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];

    // *** DIAGRAM 4724 OF 15495 ***
    // Wavefunction(s) for diagram number 4724
    // (none)
    // Amplitude(s) for diagram number 4724
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[2], w_fp[361], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[584] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[704] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4725 OF 15495 ***
    // Wavefunction(s) for diagram number 4725
    // (none)
    // Amplitude(s) for diagram number 4725
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[449], w_fp[128], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];

    // *** DIAGRAM 4726 OF 15495 ***
    // Wavefunction(s) for diagram number 4726
    // (none)
    // Amplitude(s) for diagram number 4726
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[2], w_fp[472], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[46] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[256] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[262] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4727 OF 15495 ***
    // Wavefunction(s) for diagram number 4727
    // (none)
    // Amplitude(s) for diagram number 4727
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[472], w_fp[1], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[147] -= amp_sv[0];
    jamp_sv[149] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[212] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4728 OF 15495 ***
    // Wavefunction(s) for diagram number 4728
    // (none)
    // Amplitude(s) for diagram number 4728
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[128], w_fp[435], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];

    // *** DIAGRAM 4729 OF 15495 ***
    // Wavefunction(s) for diagram number 4729
    // (none)
    // Amplitude(s) for diagram number 4729
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[361], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];

    // *** DIAGRAM 4730 OF 15495 ***
    // Wavefunction(s) for diagram number 4730
    // (none)
    // Amplitude(s) for diagram number 4730
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[2], w_fp[587], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] -= amp_sv[0];
    jamp_sv[29] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[243] += amp_sv[0];
    jamp_sv[245] -= amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[584] += amp_sv[0];
    jamp_sv[702] += amp_sv[0];
    jamp_sv[704] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[2], w_fp[586], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[2], w_fp[585], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[27] += amp_sv[0];
    jamp_sv[29] -= amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[125] += amp_sv[0];
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[243] -= amp_sv[0];
    jamp_sv[245] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[584] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[702] -= amp_sv[0];
    jamp_sv[704] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];

    // *** DIAGRAM 4731 OF 15495 ***
    // Wavefunction(s) for diagram number 4731
    // (none)
    // Amplitude(s) for diagram number 4731
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[150], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4732 OF 15495 ***
    // Wavefunction(s) for diagram number 4732
    // (none)
    // Amplitude(s) for diagram number 4732
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[7], w_fp[549], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4733 OF 15495 ***
    // Wavefunction(s) for diagram number 4733
    // (none)
    // Amplitude(s) for diagram number 4733
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[5], w_fp[550], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4734 OF 15495 ***
    // Wavefunction(s) for diagram number 4734
    // (none)
    // Amplitude(s) for diagram number 4734
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[453], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];

    // *** DIAGRAM 4735 OF 15495 ***
    // Wavefunction(s) for diagram number 4735
    // (none)
    // Amplitude(s) for diagram number 4735
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4736 OF 15495 ***
    // Wavefunction(s) for diagram number 4736
    // (none)
    // Amplitude(s) for diagram number 4736
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[453], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];

    // *** DIAGRAM 4737 OF 15495 ***
    // Wavefunction(s) for diagram number 4737
    // (none)
    // Amplitude(s) for diagram number 4737
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[549], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4738 OF 15495 ***
    // Wavefunction(s) for diagram number 4738
    // (none)
    // Amplitude(s) for diagram number 4738
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[533], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4739 OF 15495 ***
    // Wavefunction(s) for diagram number 4739
    // (none)
    // Amplitude(s) for diagram number 4739
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[532], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4740 OF 15495 ***
    // Wavefunction(s) for diagram number 4740
    // (none)
    // Amplitude(s) for diagram number 4740
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[527], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4741 OF 15495 ***
    // Wavefunction(s) for diagram number 4741
    // (none)
    // Amplitude(s) for diagram number 4741
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[532], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4742 OF 15495 ***
    // Wavefunction(s) for diagram number 4742
    // (none)
    // Amplitude(s) for diagram number 4742
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[527], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4743 OF 15495 ***
    // Wavefunction(s) for diagram number 4743
    // (none)
    // Amplitude(s) for diagram number 4743
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[533], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4744 OF 15495 ***
    // Wavefunction(s) for diagram number 4744
    // (none)
    // Amplitude(s) for diagram number 4744
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[535], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];

    // *** DIAGRAM 4745 OF 15495 ***
    // Wavefunction(s) for diagram number 4745
    // (none)
    // Amplitude(s) for diagram number 4745
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[490], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4746 OF 15495 ***
    // Wavefunction(s) for diagram number 4746
    // (none)
    // Amplitude(s) for diagram number 4746
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4747 OF 15495 ***
    // Wavefunction(s) for diagram number 4747
    // (none)
    // Amplitude(s) for diagram number 4747
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[7], w_fp[553], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4748 OF 15495 ***
    // Wavefunction(s) for diagram number 4748
    // (none)
    // Amplitude(s) for diagram number 4748
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[150], w_fp[490], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4749 OF 15495 ***
    // Wavefunction(s) for diagram number 4749
    // (none)
    // Amplitude(s) for diagram number 4749
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[553], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4750 OF 15495 ***
    // Wavefunction(s) for diagram number 4750
    // (none)
    // Amplitude(s) for diagram number 4750
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[535], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];

    // *** DIAGRAM 4751 OF 15495 ***
    // Wavefunction(s) for diagram number 4751
    // (none)
    // Amplitude(s) for diagram number 4751
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[537], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];

    // *** DIAGRAM 4752 OF 15495 ***
    // Wavefunction(s) for diagram number 4752
    // (none)
    // Amplitude(s) for diagram number 4752
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4753 OF 15495 ***
    // Wavefunction(s) for diagram number 4753
    // (none)
    // Amplitude(s) for diagram number 4753
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[437], w_fp[1], w_fp[150], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 4754 OF 15495 ***
    // Wavefunction(s) for diagram number 4754
    // (none)
    // Amplitude(s) for diagram number 4754
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[5], w_fp[575], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4755 OF 15495 ***
    // Wavefunction(s) for diagram number 4755
    // (none)
    // Amplitude(s) for diagram number 4755
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[150], w_fp[480], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];

    // *** DIAGRAM 4756 OF 15495 ***
    // Wavefunction(s) for diagram number 4756
    // (none)
    // Amplitude(s) for diagram number 4756
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4757 OF 15495 ***
    // Wavefunction(s) for diagram number 4757
    // (none)
    // Amplitude(s) for diagram number 4757
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[537], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];

    // *** DIAGRAM 4758 OF 15495 ***
    // Wavefunction(s) for diagram number 4758
    // (none)
    // Amplitude(s) for diagram number 4758
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[150], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[150], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[578], w_fp[150], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4759 OF 15495 ***
    // Wavefunction(s) for diagram number 4759
    // (none)
    // Amplitude(s) for diagram number 4759
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[577], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[578], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4760 OF 15495 ***
    // Wavefunction(s) for diagram number 4760
    // (none)
    // Amplitude(s) for diagram number 4760
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[582], w_fp[150], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += amp_sv[0];
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[150], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[150], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];

    // *** DIAGRAM 4761 OF 15495 ***
    // Wavefunction(s) for diagram number 4761
    // (none)
    // Amplitude(s) for diagram number 4761
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[582], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[583], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[584], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[606] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4762 OF 15495 ***
    // Wavefunction(s) for diagram number 4762
    // (none)
    // Amplitude(s) for diagram number 4762
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[481], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[510], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[448], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[608] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4763 OF 15495 ***
    // Wavefunction(s) for diagram number 4763
    // (none)
    // Amplitude(s) for diagram number 4763
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[510], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[448], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4764 OF 15495 ***
    // Wavefunction(s) for diagram number 4764
    // (none)
    // Amplitude(s) for diagram number 4764
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[150], w_fp[100], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4765 OF 15495 ***
    // Wavefunction(s) for diagram number 4765
    // (none)
    // Amplitude(s) for diagram number 4765
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[122], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4766 OF 15495 ***
    // Wavefunction(s) for diagram number 4766
    // (none)
    // Amplitude(s) for diagram number 4766
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[109], w_fp[2], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4767 OF 15495 ***
    // Wavefunction(s) for diagram number 4767
    // (none)
    // Amplitude(s) for diagram number 4767
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[530], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];

    // *** DIAGRAM 4768 OF 15495 ***
    // Wavefunction(s) for diagram number 4768
    // (none)
    // Amplitude(s) for diagram number 4768
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[530], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4769 OF 15495 ***
    // Wavefunction(s) for diagram number 4769
    // (none)
    // Amplitude(s) for diagram number 4769
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[109], w_fp[530], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];

    // *** DIAGRAM 4770 OF 15495 ***
    // Wavefunction(s) for diagram number 4770
    // (none)
    // Amplitude(s) for diagram number 4770
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[471], w_fp[2], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4771 OF 15495 ***
    // Wavefunction(s) for diagram number 4771
    // (none)
    // Amplitude(s) for diagram number 4771
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[471], w_fp[122], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];

    // *** DIAGRAM 4772 OF 15495 ***
    // Wavefunction(s) for diagram number 4772
    // (none)
    // Amplitude(s) for diagram number 4772
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[474], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4773 OF 15495 ***
    // Wavefunction(s) for diagram number 4773
    // (none)
    // Amplitude(s) for diagram number 4773
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[474], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4774 OF 15495 ***
    // Wavefunction(s) for diagram number 4774
    // (none)
    // Amplitude(s) for diagram number 4774
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[122], w_fp[435], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 4775 OF 15495 ***
    // Wavefunction(s) for diagram number 4775
    // (none)
    // Amplitude(s) for diagram number 4775
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[360], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 4776 OF 15495 ***
    // Wavefunction(s) for diagram number 4776
    // (none)
    // Amplitude(s) for diagram number 4776
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[521], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[520], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[519], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 4777 OF 15495 ***
    // Wavefunction(s) for diagram number 4777
    // (none)
    // Amplitude(s) for diagram number 4777
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[546], w_fp[144], w_fp[5], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 4778 OF 15495 ***
    // Wavefunction(s) for diagram number 4778
    // (none)
    // Amplitude(s) for diagram number 4778
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[6], w_fp[549], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 4779 OF 15495 ***
    // Wavefunction(s) for diagram number 4779
    // (none)
    // Amplitude(s) for diagram number 4779
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[5], w_fp[561], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[440] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 4780 OF 15495 ***
    // Wavefunction(s) for diagram number 4780
    // (none)
    // Amplitude(s) for diagram number 4780
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[453], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[242] -= amp_sv[0];

    // *** DIAGRAM 4781 OF 15495 ***
    // Wavefunction(s) for diagram number 4781
    // (none)
    // Amplitude(s) for diagram number 4781
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[2], w_fp[561], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[486] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[510] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4782 OF 15495 ***
    // Wavefunction(s) for diagram number 4782
    // (none)
    // Amplitude(s) for diagram number 4782
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[453], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];

    // *** DIAGRAM 4783 OF 15495 ***
    // Wavefunction(s) for diagram number 4783
    // (none)
    // Amplitude(s) for diagram number 4783
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[549], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4784 OF 15495 ***
    // Wavefunction(s) for diagram number 4784
    // (none)
    // Amplitude(s) for diagram number 4784
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[533], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4785 OF 15495 ***
    // Wavefunction(s) for diagram number 4785
    // (none)
    // Amplitude(s) for diagram number 4785
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[531], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4786 OF 15495 ***
    // Wavefunction(s) for diagram number 4786
    // (none)
    // Amplitude(s) for diagram number 4786
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[527], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[26] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4787 OF 15495 ***
    // Wavefunction(s) for diagram number 4787
    // (none)
    // Amplitude(s) for diagram number 4787
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[531], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4788 OF 15495 ***
    // Wavefunction(s) for diagram number 4788
    // (none)
    // Amplitude(s) for diagram number 4788
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[527], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4789 OF 15495 ***
    // Wavefunction(s) for diagram number 4789
    // (none)
    // Amplitude(s) for diagram number 4789
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[533], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4790 OF 15495 ***
    // Wavefunction(s) for diagram number 4790
    // (none)
    // Amplitude(s) for diagram number 4790
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[535], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[248] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];

    // *** DIAGRAM 4791 OF 15495 ***
    // Wavefunction(s) for diagram number 4791
    // (none)
    // Amplitude(s) for diagram number 4791
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[248] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[554] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4792 OF 15495 ***
    // Wavefunction(s) for diagram number 4792
    // (none)
    // Amplitude(s) for diagram number 4792
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[451], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 4793 OF 15495 ***
    // Wavefunction(s) for diagram number 4793
    // (none)
    // Amplitude(s) for diagram number 4793
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[6], w_fp[553], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[144] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[510] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[554] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 4794 OF 15495 ***
    // Wavefunction(s) for diagram number 4794
    // (none)
    // Amplitude(s) for diagram number 4794
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[144], w_fp[514], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[144] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[204] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[248] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[554] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 4795 OF 15495 ***
    // Wavefunction(s) for diagram number 4795
    // (none)
    // Amplitude(s) for diagram number 4795
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[553], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4796 OF 15495 ***
    // Wavefunction(s) for diagram number 4796
    // (none)
    // Amplitude(s) for diagram number 4796
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[535], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];

    // *** DIAGRAM 4797 OF 15495 ***
    // Wavefunction(s) for diagram number 4797
    // (none)
    // Amplitude(s) for diagram number 4797
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[536], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];

    // *** DIAGRAM 4798 OF 15495 ***
    // Wavefunction(s) for diagram number 4798
    // (none)
    // Amplitude(s) for diagram number 4798
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[497], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[254] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[434] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 4799 OF 15495 ***
    // Wavefunction(s) for diagram number 4799
    // (none)
    // Amplitude(s) for diagram number 4799
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[488], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 4800 OF 15495 ***
    // Wavefunction(s) for diagram number 4800
    // (none)
    // Amplitude(s) for diagram number 4800
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[5], w_fp[571], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[146] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[434] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

#ifdef MGONGPUCPP_GPUIMPL
    // *** STORE JAMPS ***
    // In CUDA, copy the local jamp to the output global-memory jamp
    constexpr int ihel0 = 0; // the allJamps buffer already points to a specific helicity _within a super-buffer for nGoodHel helicities_
    for( int icol = 0; icol < ncolor; icol++ )
      J_ACCESS::kernelAccessIcolIhelNhel( jamps, icol, ihel0, nGoodHel ) += jamp_sv[icol]; // update jamps
#else
    // In C++, copy the local jamp to the output array passed as function argument
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
