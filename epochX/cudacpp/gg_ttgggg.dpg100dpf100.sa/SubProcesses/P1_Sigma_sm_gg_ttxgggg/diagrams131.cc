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
  diagramgroup131( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 99 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 123 );
    retrieveWf( wfs, w_cx, nevt, 126 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 143 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 171 );
    retrieveWf( wfs, w_cx, nevt, 178 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 197 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 211 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 217 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 219 );
    retrieveWf( wfs, w_cx, nevt, 220 );
    retrieveWf( wfs, w_cx, nevt, 222 );
    retrieveWf( wfs, w_cx, nevt, 223 );
    retrieveWf( wfs, w_cx, nevt, 224 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 256 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 450 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 489 );
    retrieveWf( wfs, w_cx, nevt, 509 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 658 );
    retrieveWf( wfs, w_cx, nevt, 660 );
    retrieveWf( wfs, w_cx, nevt, 674 );
    retrieveWf( wfs, w_cx, nevt, 684 );
    retrieveWf( wfs, w_cx, nevt, 685 );
    retrieveWf( wfs, w_cx, nevt, 686 );
    retrieveWf( wfs, w_cx, nevt, 687 );
    retrieveWf( wfs, w_cx, nevt, 688 );
    retrieveWf( wfs, w_cx, nevt, 689 );
#endif
#endif

    // *** DIAGRAM 13001 OF 15495 ***
    // Wavefunction(s) for diagram number 13001
    // (none)
    // Amplitude(s) for diagram number 13001
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[523] -= amp_sv[0];

    // *** DIAGRAM 13002 OF 15495 ***
    // Wavefunction(s) for diagram number 13002
    // (none)
    // Amplitude(s) for diagram number 13002
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[40], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[507] -= amp_sv[0];

    // *** DIAGRAM 13003 OF 15495 ***
    // Wavefunction(s) for diagram number 13003
    // (none)
    // Amplitude(s) for diagram number 13003
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[517] -= amp_sv[0];

    // *** DIAGRAM 13004 OF 15495 ***
    // Wavefunction(s) for diagram number 13004
    // (none)
    // Amplitude(s) for diagram number 13004
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[674], w_fp[212], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] -= amp_sv[0];

    // *** DIAGRAM 13005 OF 15495 ***
    // Wavefunction(s) for diagram number 13005
    // (none)
    // Amplitude(s) for diagram number 13005
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[687], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[556] -= amp_sv[0];

    // *** DIAGRAM 13006 OF 15495 ***
    // Wavefunction(s) for diagram number 13006
    // (none)
    // Amplitude(s) for diagram number 13006
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[674], w_fp[213], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[594] -= amp_sv[0];

    // *** DIAGRAM 13007 OF 15495 ***
    // Wavefunction(s) for diagram number 13007
    // (none)
    // Amplitude(s) for diagram number 13007
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[688], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] -= amp_sv[0];

    // *** DIAGRAM 13008 OF 15495 ***
    // Wavefunction(s) for diagram number 13008
    // (none)
    // Amplitude(s) for diagram number 13008
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[199], w_fp[687], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[553] -= amp_sv[0];

    // *** DIAGRAM 13009 OF 15495 ***
    // Wavefunction(s) for diagram number 13009
    // (none)
    // Amplitude(s) for diagram number 13009
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[532], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[559] -= amp_sv[0];

    // *** DIAGRAM 13010 OF 15495 ***
    // Wavefunction(s) for diagram number 13010
    // (none)
    // Amplitude(s) for diagram number 13010
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[216], w_fp[688], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= amp_sv[0];

    // *** DIAGRAM 13011 OF 15495 ***
    // Wavefunction(s) for diagram number 13011
    // (none)
    // Amplitude(s) for diagram number 13011
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[588], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[583] -= amp_sv[0];

    // *** DIAGRAM 13012 OF 15495 ***
    // Wavefunction(s) for diagram number 13012
    // (none)
    // Amplitude(s) for diagram number 13012
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[540], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13013 OF 15495 ***
    // Wavefunction(s) for diagram number 13013
    // (none)
    // Amplitude(s) for diagram number 13013
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[540], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[502] += amp_sv[0];

    // *** DIAGRAM 13014 OF 15495 ***
    // Wavefunction(s) for diagram number 13014
    // (none)
    // Amplitude(s) for diagram number 13014
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[123], w_fp[540], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13015 OF 15495 ***
    // Wavefunction(s) for diagram number 13015
    // (none)
    // Amplitude(s) for diagram number 13015
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[476], w_fp[100], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13016 OF 15495 ***
    // Wavefunction(s) for diagram number 13016
    // (none)
    // Amplitude(s) for diagram number 13016
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[197], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[520] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[572] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];

    // *** DIAGRAM 13017 OF 15495 ***
    // Wavefunction(s) for diagram number 13017
    // (none)
    // Amplitude(s) for diagram number 13017
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[224], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13018 OF 15495 ***
    // Wavefunction(s) for diagram number 13018
    // (none)
    // Amplitude(s) for diagram number 13018
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[476], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[507] += amp_sv[0];
    jamp_sv[509] -= amp_sv[0];
    jamp_sv[520] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];

    // *** DIAGRAM 13019 OF 15495 ***
    // Wavefunction(s) for diagram number 13019
    // (none)
    // Amplitude(s) for diagram number 13019
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[197], w_fp[522], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += amp_sv[0];
    jamp_sv[502] -= amp_sv[0];
    jamp_sv[570] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];

    // *** DIAGRAM 13020 OF 15495 ***
    // Wavefunction(s) for diagram number 13020
    // (none)
    // Amplitude(s) for diagram number 13020
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[217], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13021 OF 15495 ***
    // Wavefunction(s) for diagram number 13021
    // (none)
    // Amplitude(s) for diagram number 13021
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[123], w_fp[476], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13022 OF 15495 ***
    // Wavefunction(s) for diagram number 13022
    // (none)
    // Amplitude(s) for diagram number 13022
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[355], w_fp[224], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[570] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13023 OF 15495 ***
    // Wavefunction(s) for diagram number 13023
    // (none)
    // Amplitude(s) for diagram number 13023
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[217], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13024 OF 15495 ***
    // Wavefunction(s) for diagram number 13024
    // (none)
    // Amplitude(s) for diagram number 13024
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[197], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[197], w_fp[126], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[502] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[520] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[526] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[197], w_fp[62], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[485] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[570] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[572] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13025 OF 15495 ***
    // Wavefunction(s) for diagram number 13025
    // (none)
    // Amplitude(s) for diagram number 13025
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[686], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] -= amp_sv[0];

    // *** DIAGRAM 13026 OF 15495 ***
    // Wavefunction(s) for diagram number 13026
    // (none)
    // Amplitude(s) for diagram number 13026
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[685], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[500] -= amp_sv[0];

    // *** DIAGRAM 13027 OF 15495 ***
    // Wavefunction(s) for diagram number 13027
    // (none)
    // Amplitude(s) for diagram number 13027
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[12], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[484] -= amp_sv[0];

    // *** DIAGRAM 13028 OF 15495 ***
    // Wavefunction(s) for diagram number 13028
    // (none)
    // Amplitude(s) for diagram number 13028
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[685], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] -= amp_sv[0];

    // *** DIAGRAM 13029 OF 15495 ***
    // Wavefunction(s) for diagram number 13029
    // (none)
    // Amplitude(s) for diagram number 13029
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[12], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= amp_sv[0];

    // *** DIAGRAM 13030 OF 15495 ***
    // Wavefunction(s) for diagram number 13030
    // (none)
    // Amplitude(s) for diagram number 13030
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[686], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[487] -= amp_sv[0];

    // *** DIAGRAM 13031 OF 15495 ***
    // Wavefunction(s) for diagram number 13031
    // (none)
    // Amplitude(s) for diagram number 13031
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[524] -= amp_sv[0];

    // *** DIAGRAM 13032 OF 15495 ***
    // Wavefunction(s) for diagram number 13032
    // (none)
    // Amplitude(s) for diagram number 13032
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[543], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[514] -= amp_sv[0];

    // *** DIAGRAM 13033 OF 15495 ***
    // Wavefunction(s) for diagram number 13033
    // (none)
    // Amplitude(s) for diagram number 13033
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[211], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[548] -= amp_sv[0];

    // *** DIAGRAM 13034 OF 15495 ***
    // Wavefunction(s) for diagram number 13034
    // (none)
    // Amplitude(s) for diagram number 13034
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[543], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[538] -= amp_sv[0];

    // *** DIAGRAM 13035 OF 15495 ***
    // Wavefunction(s) for diagram number 13035
    // (none)
    // Amplitude(s) for diagram number 13035
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[660], w_fp[213], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[590] -= amp_sv[0];

    // *** DIAGRAM 13036 OF 15495 ***
    // Wavefunction(s) for diagram number 13036
    // (none)
    // Amplitude(s) for diagram number 13036
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[509], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[584] -= amp_sv[0];

    // *** DIAGRAM 13037 OF 15495 ***
    // Wavefunction(s) for diagram number 13037
    // (none)
    // Amplitude(s) for diagram number 13037
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[40], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[508] -= amp_sv[0];

    // *** DIAGRAM 13038 OF 15495 ***
    // Wavefunction(s) for diagram number 13038
    // (none)
    // Amplitude(s) for diagram number 13038
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[476], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[522] -= amp_sv[0];

    // *** DIAGRAM 13039 OF 15495 ***
    // Wavefunction(s) for diagram number 13039
    // (none)
    // Amplitude(s) for diagram number 13039
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[40], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] -= amp_sv[0];

    // *** DIAGRAM 13040 OF 15495 ***
    // Wavefunction(s) for diagram number 13040
    // (none)
    // Amplitude(s) for diagram number 13040
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[511] -= amp_sv[0];

    // *** DIAGRAM 13041 OF 15495 ***
    // Wavefunction(s) for diagram number 13041
    // (none)
    // Amplitude(s) for diagram number 13041
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[658], w_fp[211], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[546] -= amp_sv[0];

    // *** DIAGRAM 13042 OF 15495 ***
    // Wavefunction(s) for diagram number 13042
    // (none)
    // Amplitude(s) for diagram number 13042
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[689], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[532] -= amp_sv[0];

    // *** DIAGRAM 13043 OF 15495 ***
    // Wavefunction(s) for diagram number 13043
    // (none)
    // Amplitude(s) for diagram number 13043
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[658], w_fp[213], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] -= amp_sv[0];

    // *** DIAGRAM 13044 OF 15495 ***
    // Wavefunction(s) for diagram number 13044
    // (none)
    // Amplitude(s) for diagram number 13044
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[688], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] -= amp_sv[0];

    // *** DIAGRAM 13045 OF 15495 ***
    // Wavefunction(s) for diagram number 13045
    // (none)
    // Amplitude(s) for diagram number 13045
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[171], w_fp[689], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[529] -= amp_sv[0];

    // *** DIAGRAM 13046 OF 15495 ***
    // Wavefunction(s) for diagram number 13046
    // (none)
    // Amplitude(s) for diagram number 13046
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[450], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[535] -= amp_sv[0];

    // *** DIAGRAM 13047 OF 15495 ***
    // Wavefunction(s) for diagram number 13047
    // (none)
    // Amplitude(s) for diagram number 13047
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[218], w_fp[688], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[576] -= amp_sv[0];

    // *** DIAGRAM 13048 OF 15495 ***
    // Wavefunction(s) for diagram number 13048
    // (none)
    // Amplitude(s) for diagram number 13048
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[481], w_fp[213], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[582] -= amp_sv[0];

    // *** DIAGRAM 13049 OF 15495 ***
    // Wavefunction(s) for diagram number 13049
    // (none)
    // Amplitude(s) for diagram number 13049
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[540], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13050 OF 15495 ***
    // Wavefunction(s) for diagram number 13050
    // (none)
    // Amplitude(s) for diagram number 13050
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[540], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[500] += amp_sv[0];

    // *** DIAGRAM 13051 OF 15495 ***
    // Wavefunction(s) for diagram number 13051
    // (none)
    // Amplitude(s) for diagram number 13051
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[540], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13052 OF 15495 ***
    // Wavefunction(s) for diagram number 13052
    // (none)
    // Amplitude(s) for diagram number 13052
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[476], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13053 OF 15495 ***
    // Wavefunction(s) for diagram number 13053
    // (none)
    // Amplitude(s) for diagram number 13053
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[197], w_fp[363], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[514] += amp_sv[0];
    jamp_sv[524] -= amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[590] += amp_sv[0];

    // *** DIAGRAM 13054 OF 15495 ***
    // Wavefunction(s) for diagram number 13054
    // (none)
    // Amplitude(s) for diagram number 13054
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[223], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13055 OF 15495 ***
    // Wavefunction(s) for diagram number 13055
    // (none)
    // Amplitude(s) for diagram number 13055
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[476], w_fp[518], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] += amp_sv[0];
    jamp_sv[508] -= amp_sv[0];
    jamp_sv[514] -= amp_sv[0];
    jamp_sv[524] += amp_sv[0];

    // *** DIAGRAM 13056 OF 15495 ***
    // Wavefunction(s) for diagram number 13056
    // (none)
    // Amplitude(s) for diagram number 13056
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[197], w_fp[518], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] += amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];

    // *** DIAGRAM 13057 OF 15495 ***
    // Wavefunction(s) for diagram number 13057
    // (none)
    // Amplitude(s) for diagram number 13057
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[1], w_fp[219], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13058 OF 15495 ***
    // Wavefunction(s) for diagram number 13058
    // (none)
    // Amplitude(s) for diagram number 13058
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[99], w_fp[476], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13059 OF 15495 ***
    // Wavefunction(s) for diagram number 13059
    // (none)
    // Amplitude(s) for diagram number 13059
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[256], w_fp[223], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13060 OF 15495 ***
    // Wavefunction(s) for diagram number 13060
    // (none)
    // Amplitude(s) for diagram number 13060
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[363], w_fp[219], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13061 OF 15495 ***
    // Wavefunction(s) for diagram number 13061
    // (none)
    // Amplitude(s) for diagram number 13061
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[19], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[18], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[490] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[500] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[514] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[524] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[197], w_fp[17], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[481] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[484] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[590] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13062 OF 15495 ***
    // Wavefunction(s) for diagram number 13062
    // (none)
    // Amplitude(s) for diagram number 13062
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[686], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= amp_sv[0];

    // *** DIAGRAM 13063 OF 15495 ***
    // Wavefunction(s) for diagram number 13063
    // (none)
    // Amplitude(s) for diagram number 13063
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[684], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[494] -= amp_sv[0];

    // *** DIAGRAM 13064 OF 15495 ***
    // Wavefunction(s) for diagram number 13064
    // (none)
    // Amplitude(s) for diagram number 13064
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[12], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[482] -= amp_sv[0];

    // *** DIAGRAM 13065 OF 15495 ***
    // Wavefunction(s) for diagram number 13065
    // (none)
    // Amplitude(s) for diagram number 13065
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[684], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[492] -= amp_sv[0];

    // *** DIAGRAM 13066 OF 15495 ***
    // Wavefunction(s) for diagram number 13066
    // (none)
    // Amplitude(s) for diagram number 13066
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[12], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= amp_sv[0];

    // *** DIAGRAM 13067 OF 15495 ***
    // Wavefunction(s) for diagram number 13067
    // (none)
    // Amplitude(s) for diagram number 13067
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[686], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[486] -= amp_sv[0];

    // *** DIAGRAM 13068 OF 15495 ***
    // Wavefunction(s) for diagram number 13068
    // (none)
    // Amplitude(s) for diagram number 13068
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[518] -= amp_sv[0];

    // *** DIAGRAM 13069 OF 15495 ***
    // Wavefunction(s) for diagram number 13069
    // (none)
    // Amplitude(s) for diagram number 13069
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] -= amp_sv[0];

    // *** DIAGRAM 13070 OF 15495 ***
    // Wavefunction(s) for diagram number 13070
    // (none)
    // Amplitude(s) for diagram number 13070
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[178], w_fp[211], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 13071 OF 15495 ***
    // Wavefunction(s) for diagram number 13071
    // (none)
    // Amplitude(s) for diagram number 13071
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[536] -= amp_sv[0];

    // *** DIAGRAM 13072 OF 15495 ***
    // Wavefunction(s) for diagram number 13072
    // (none)
    // Amplitude(s) for diagram number 13072
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[178], w_fp[212], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[566] -= amp_sv[0];

    // *** DIAGRAM 13073 OF 15495 ***
    // Wavefunction(s) for diagram number 13073
    // (none)
    // Amplitude(s) for diagram number 13073
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 13074 OF 15495 ***
    // Wavefunction(s) for diagram number 13074
    // (none)
    // Amplitude(s) for diagram number 13074
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[40], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[506] -= amp_sv[0];

    // *** DIAGRAM 13075 OF 15495 ***
    // Wavefunction(s) for diagram number 13075
    // (none)
    // Amplitude(s) for diagram number 13075
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[476], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[516] -= amp_sv[0];

    // *** DIAGRAM 13076 OF 15495 ***
    // Wavefunction(s) for diagram number 13076
    // (none)
    // Amplitude(s) for diagram number 13076
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[40], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] -= amp_sv[0];

    // *** DIAGRAM 13077 OF 15495 ***
    // Wavefunction(s) for diagram number 13077
    // (none)
    // Amplitude(s) for diagram number 13077
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[476], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[510] -= amp_sv[0];

    // *** DIAGRAM 13078 OF 15495 ***
    // Wavefunction(s) for diagram number 13078
    // (none)
    // Amplitude(s) for diagram number 13078
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[187], w_fp[211], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13079 OF 15495 ***
    // Wavefunction(s) for diagram number 13079
    // (none)
    // Amplitude(s) for diagram number 13079
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[689], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 13080 OF 15495 ***
    // Wavefunction(s) for diagram number 13080
    // (none)
    // Amplitude(s) for diagram number 13080
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[187], w_fp[212], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 13081 OF 15495 ***
    // Wavefunction(s) for diagram number 13081
    // (none)
    // Amplitude(s) for diagram number 13081
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[687], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[554] -= amp_sv[0];

    // *** DIAGRAM 13082 OF 15495 ***
    // Wavefunction(s) for diagram number 13082
    // (none)
    // Amplitude(s) for diagram number 13082
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[689], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[528] -= amp_sv[0];

    // *** DIAGRAM 13083 OF 15495 ***
    // Wavefunction(s) for diagram number 13083
    // (none)
    // Amplitude(s) for diagram number 13083
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[211], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 13084 OF 15495 ***
    // Wavefunction(s) for diagram number 13084
    // (none)
    // Amplitude(s) for diagram number 13084
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[687], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[552] -= amp_sv[0];

    // *** DIAGRAM 13085 OF 15495 ***
    // Wavefunction(s) for diagram number 13085
    // (none)
    // Amplitude(s) for diagram number 13085
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[212], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 13086 OF 15495 ***
    // Wavefunction(s) for diagram number 13086
    // (none)
    // Amplitude(s) for diagram number 13086
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[540], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13087 OF 15495 ***
    // Wavefunction(s) for diagram number 13087
    // (none)
    // Amplitude(s) for diagram number 13087
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[540], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];

    // *** DIAGRAM 13088 OF 15495 ***
    // Wavefunction(s) for diagram number 13088
    // (none)
    // Amplitude(s) for diagram number 13088
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[143], w_fp[540], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13089 OF 15495 ***
    // Wavefunction(s) for diagram number 13089
    // (none)
    // Amplitude(s) for diagram number 13089
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[476], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13090 OF 15495 ***
    // Wavefunction(s) for diagram number 13090
    // (none)
    // Amplitude(s) for diagram number 13090
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[197], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13091 OF 15495 ***
    // Wavefunction(s) for diagram number 13091
    // (none)
    // Amplitude(s) for diagram number 13091
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[222], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13092 OF 15495 ***
    // Wavefunction(s) for diagram number 13092
    // (none)
    // Amplitude(s) for diagram number 13092
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[476], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];

    // *** DIAGRAM 13093 OF 15495 ***
    // Wavefunction(s) for diagram number 13093
    // (none)
    // Amplitude(s) for diagram number 13093
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[197], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13094 OF 15495 ***
    // Wavefunction(s) for diagram number 13094
    // (none)
    // Amplitude(s) for diagram number 13094
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[220], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13095 OF 15495 ***
    // Wavefunction(s) for diagram number 13095
    // (none)
    // Amplitude(s) for diagram number 13095
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[143], w_fp[476], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13096 OF 15495 ***
    // Wavefunction(s) for diagram number 13096
    // (none)
    // Amplitude(s) for diagram number 13096
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[222], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13097 OF 15495 ***
    // Wavefunction(s) for diagram number 13097
    // (none)
    // Amplitude(s) for diagram number 13097
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[220], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13098 OF 15495 ***
    // Wavefunction(s) for diagram number 13098
    // (none)
    // Amplitude(s) for diagram number 13098
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[197], w_fp[136], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[197], w_fp[137], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[197], w_fp[551], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13099 OF 15495 ***
    // Wavefunction(s) for diagram number 13099
    // (none)
    // Amplitude(s) for diagram number 13099
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[254], w_fp[7], w_fp[489], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[480] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[498] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[499] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[501] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[503] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13100 OF 15495 ***
    // Wavefunction(s) for diagram number 13100
    // (none)
    // Amplitude(s) for diagram number 13100
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[685], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[503] += amp_sv[0];

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
