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
  diagramgroup120( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 0 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 149 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 154 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 200 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 317 );
    retrieveWf( wfs, w_cx, nevt, 325 );
    retrieveWf( wfs, w_cx, nevt, 328 );
    retrieveWf( wfs, w_cx, nevt, 329 );
    retrieveWf( wfs, w_cx, nevt, 332 );
    retrieveWf( wfs, w_cx, nevt, 333 );
    retrieveWf( wfs, w_cx, nevt, 334 );
    retrieveWf( wfs, w_cx, nevt, 335 );
    retrieveWf( wfs, w_cx, nevt, 336 );
    retrieveWf( wfs, w_cx, nevt, 337 );
    retrieveWf( wfs, w_cx, nevt, 344 );
    retrieveWf( wfs, w_cx, nevt, 347 );
    retrieveWf( wfs, w_cx, nevt, 348 );
    retrieveWf( wfs, w_cx, nevt, 349 );
    retrieveWf( wfs, w_cx, nevt, 350 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 501 );
    retrieveWf( wfs, w_cx, nevt, 503 );
    retrieveWf( wfs, w_cx, nevt, 504 );
    retrieveWf( wfs, w_cx, nevt, 507 );
    retrieveWf( wfs, w_cx, nevt, 508 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 548 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 554 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 597 );
    retrieveWf( wfs, w_cx, nevt, 599 );
    retrieveWf( wfs, w_cx, nevt, 601 );
    retrieveWf( wfs, w_cx, nevt, 603 );
    retrieveWf( wfs, w_cx, nevt, 605 );
    retrieveWf( wfs, w_cx, nevt, 607 );
    retrieveWf( wfs, w_cx, nevt, 608 );
    retrieveWf( wfs, w_cx, nevt, 610 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 696 );
    retrieveWf( wfs, w_cx, nevt, 698 );
    retrieveWf( wfs, w_cx, nevt, 703 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 707 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 711 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 724 );
    retrieveWf( wfs, w_cx, nevt, 725 );
    retrieveWf( wfs, w_cx, nevt, 743 );
#endif
#endif

    // *** DIAGRAM 11901 OF 15495 ***
    // Wavefunction(s) for diagram number 11901
    // (none)
    // Amplitude(s) for diagram number 11901
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[501], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

    // *** DIAGRAM 11902 OF 15495 ***
    // Wavefunction(s) for diagram number 11902
    // (none)
    // Amplitude(s) for diagram number 11902
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[2], w_fp[350], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[226] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[322] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[332] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[538] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[548] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11903 OF 15495 ***
    // Wavefunction(s) for diagram number 11903
    // (none)
    // Amplitude(s) for diagram number 11903
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[527], w_fp[118], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];

    // *** DIAGRAM 11904 OF 15495 ***
    // Wavefunction(s) for diagram number 11904
    // (none)
    // Amplitude(s) for diagram number 11904
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[501], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11905 OF 15495 ***
    // Wavefunction(s) for diagram number 11905
    // (none)
    // Amplitude(s) for diagram number 11905
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[347], w_fp[2], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[40] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[79] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[313] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[529] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11906 OF 15495 ***
    // Wavefunction(s) for diagram number 11906
    // (none)
    // Amplitude(s) for diagram number 11906
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[325], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[226] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[644] -= amp_sv[0];

    // *** DIAGRAM 11907 OF 15495 ***
    // Wavefunction(s) for diagram number 11907
    // (none)
    // Amplitude(s) for diagram number 11907
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[200], w_fp[501], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 11908 OF 15495 ***
    // Wavefunction(s) for diagram number 11908
    // (none)
    // Amplitude(s) for diagram number 11908
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[347], w_fp[118], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];

    // *** DIAGRAM 11909 OF 15495 ***
    // Wavefunction(s) for diagram number 11909
    // (none)
    // Amplitude(s) for diagram number 11909
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[350], w_fp[154], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

    // *** DIAGRAM 11910 OF 15495 ***
    // Wavefunction(s) for diagram number 11910
    // (none)
    // Amplitude(s) for diagram number 11910
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[2], w_fp[707], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[322] -= amp_sv[0];
    jamp_sv[332] += amp_sv[0];
    jamp_sv[538] += amp_sv[0];
    jamp_sv[548] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[2], w_fp[329], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] += amp_sv[0];
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[79] -= amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[168], w_fp[2], w_fp[328], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] += amp_sv[0];
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[97] -= amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];

    // *** DIAGRAM 11911 OF 15495 ***
    // Wavefunction(s) for diagram number 11911
    // (none)
    // Amplitude(s) for diagram number 11911
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 11912 OF 15495 ***
    // Wavefunction(s) for diagram number 11912
    // (none)
    // Amplitude(s) for diagram number 11912
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[5], w_fp[548], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];

    // *** DIAGRAM 11913 OF 15495 ***
    // Wavefunction(s) for diagram number 11913
    // (none)
    // Amplitude(s) for diagram number 11913
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[4], w_fp[706], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 11914 OF 15495 ***
    // Wavefunction(s) for diagram number 11914
    // (none)
    // Amplitude(s) for diagram number 11914
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[508], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 11915 OF 15495 ***
    // Wavefunction(s) for diagram number 11915
    // (none)
    // Amplitude(s) for diagram number 11915
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11916 OF 15495 ***
    // Wavefunction(s) for diagram number 11916
    // (none)
    // Amplitude(s) for diagram number 11916
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[508], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];

    // *** DIAGRAM 11917 OF 15495 ***
    // Wavefunction(s) for diagram number 11917
    // (none)
    // Amplitude(s) for diagram number 11917
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[548], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11918 OF 15495 ***
    // Wavefunction(s) for diagram number 11918
    // (none)
    // Amplitude(s) for diagram number 11918
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[595], w_fp[501], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11919 OF 15495 ***
    // Wavefunction(s) for diagram number 11919
    // (none)
    // Amplitude(s) for diagram number 11919
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[501], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11920 OF 15495 ***
    // Wavefunction(s) for diagram number 11920
    // (none)
    // Amplitude(s) for diagram number 11920
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[332], w_fp[5], w_fp[724], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11921 OF 15495 ***
    // Wavefunction(s) for diagram number 11921
    // (none)
    // Amplitude(s) for diagram number 11921
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[577], w_fp[2], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];

    // *** DIAGRAM 11922 OF 15495 ***
    // Wavefunction(s) for diagram number 11922
    // (none)
    // Amplitude(s) for diagram number 11922
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[333], w_fp[4], w_fp[724], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11923 OF 15495 ***
    // Wavefunction(s) for diagram number 11923
    // (none)
    // Amplitude(s) for diagram number 11923
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[595], w_fp[2], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[230] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];

    // *** DIAGRAM 11924 OF 15495 ***
    // Wavefunction(s) for diagram number 11924
    // (none)
    // Amplitude(s) for diagram number 11924
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[335], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[336], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[337], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11925 OF 15495 ***
    // Wavefunction(s) for diagram number 11925
    // (none)
    // Amplitude(s) for diagram number 11925
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[507], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11926 OF 15495 ***
    // Wavefunction(s) for diagram number 11926
    // (none)
    // Amplitude(s) for diagram number 11926
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[501], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11927 OF 15495 ***
    // Wavefunction(s) for diagram number 11927
    // (none)
    // Amplitude(s) for diagram number 11927
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[507], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11928 OF 15495 ***
    // Wavefunction(s) for diagram number 11928
    // (none)
    // Amplitude(s) for diagram number 11928
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[501], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11929 OF 15495 ***
    // Wavefunction(s) for diagram number 11929
    // (none)
    // Amplitude(s) for diagram number 11929
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[150], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[150], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[332], w_fp[150], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];

    // *** DIAGRAM 11930 OF 15495 ***
    // Wavefunction(s) for diagram number 11930
    // (none)
    // Amplitude(s) for diagram number 11930
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[5], w_fp[743], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];

    // *** DIAGRAM 11931 OF 15495 ***
    // Wavefunction(s) for diagram number 11931
    // (none)
    // Amplitude(s) for diagram number 11931
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[332], w_fp[5], w_fp[725], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];

    // *** DIAGRAM 11932 OF 15495 ***
    // Wavefunction(s) for diagram number 11932
    // (none)
    // Amplitude(s) for diagram number 11932
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[743], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11933 OF 15495 ***
    // Wavefunction(s) for diagram number 11933
    // (none)
    // Amplitude(s) for diagram number 11933
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[530], w_fp[2], w_fp[332], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[222] += amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];

    // *** DIAGRAM 11934 OF 15495 ***
    // Wavefunction(s) for diagram number 11934
    // (none)
    // Amplitude(s) for diagram number 11934
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[333], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[333], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
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
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[333], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 11935 OF 15495 ***
    // Wavefunction(s) for diagram number 11935
    // (none)
    // Amplitude(s) for diagram number 11935
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[4], w_fp[270], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];

    // *** DIAGRAM 11936 OF 15495 ***
    // Wavefunction(s) for diagram number 11936
    // (none)
    // Amplitude(s) for diagram number 11936
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[333], w_fp[4], w_fp[725], COUPs[0], 1.0, &amp_fp[0] );
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
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 11937 OF 15495 ***
    // Wavefunction(s) for diagram number 11937
    // (none)
    // Amplitude(s) for diagram number 11937
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[270], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11938 OF 15495 ***
    // Wavefunction(s) for diagram number 11938
    // (none)
    // Amplitude(s) for diagram number 11938
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[2], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];

    // *** DIAGRAM 11939 OF 15495 ***
    // Wavefunction(s) for diagram number 11939
    // (none)
    // Amplitude(s) for diagram number 11939
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[551], w_fp[150], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[504], w_fp[150], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[222] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[544], w_fp[150], w_fp[5], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];

    // *** DIAGRAM 11940 OF 15495 ***
    // Wavefunction(s) for diagram number 11940
    // (none)
    // Amplitude(s) for diagram number 11940
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[551], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[504], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[28] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[175], w_fp[2], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[244] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[258] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[282] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[342] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11941 OF 15495 ***
    // Wavefunction(s) for diagram number 11941
    // (none)
    // Amplitude(s) for diagram number 11941
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[709], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[703], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];

    // *** DIAGRAM 11942 OF 15495 ***
    // Wavefunction(s) for diagram number 11942
    // (none)
    // Amplitude(s) for diagram number 11942
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[703], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[670], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11943 OF 15495 ***
    // Wavefunction(s) for diagram number 11943
    // (none)
    // Amplitude(s) for diagram number 11943
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[335], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[336], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
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
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[337], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[28] -= amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[66] += amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];

    // *** DIAGRAM 11944 OF 15495 ***
    // Wavefunction(s) for diagram number 11944
    // (none)
    // Amplitude(s) for diagram number 11944
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[66], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];

    // *** DIAGRAM 11945 OF 15495 ***
    // Wavefunction(s) for diagram number 11945
    // (none)
    // Amplitude(s) for diagram number 11945
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[148], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11946 OF 15495 ***
    // Wavefunction(s) for diagram number 11946
    // (none)
    // Amplitude(s) for diagram number 11946
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[149], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[96] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11947 OF 15495 ***
    // Wavefunction(s) for diagram number 11947
    // (none)
    // Amplitude(s) for diagram number 11947
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[501], w_fp[66], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 11948 OF 15495 ***
    // Wavefunction(s) for diagram number 11948
    // (none)
    // Amplitude(s) for diagram number 11948
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[349], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[418] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11949 OF 15495 ***
    // Wavefunction(s) for diagram number 11949
    // (none)
    // Amplitude(s) for diagram number 11949
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[148], w_fp[325], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];

    // *** DIAGRAM 11950 OF 15495 ***
    // Wavefunction(s) for diagram number 11950
    // (none)
    // Amplitude(s) for diagram number 11950
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[501], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11951 OF 15495 ***
    // Wavefunction(s) for diagram number 11951
    // (none)
    // Amplitude(s) for diagram number 11951
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[348], w_fp[2], w_fp[514], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11952 OF 15495 ***
    // Wavefunction(s) for diagram number 11952
    // (none)
    // Amplitude(s) for diagram number 11952
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[325], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];

    // *** DIAGRAM 11953 OF 15495 ***
    // Wavefunction(s) for diagram number 11953
    // (none)
    // Amplitude(s) for diagram number 11953
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[149], w_fp[501], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];

    // *** DIAGRAM 11954 OF 15495 ***
    // Wavefunction(s) for diagram number 11954
    // (none)
    // Amplitude(s) for diagram number 11954
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[348], w_fp[148], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];

    // *** DIAGRAM 11955 OF 15495 ***
    // Wavefunction(s) for diagram number 11955
    // (none)
    // Amplitude(s) for diagram number 11955
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[349], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];

    // *** DIAGRAM 11956 OF 15495 ***
    // Wavefunction(s) for diagram number 11956
    // (none)
    // Amplitude(s) for diagram number 11956
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[268], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[267], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[96] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];

    // *** DIAGRAM 11957 OF 15495 ***
    // Wavefunction(s) for diagram number 11957
    // (none)
    // Amplitude(s) for diagram number 11957
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[554], w_fp[148], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[303] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[305] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[423] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11958 OF 15495 ***
    // Wavefunction(s) for diagram number 11958
    // (none)
    // Amplitude(s) for diagram number 11958
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[148], w_fp[500], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[418] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];

    // *** DIAGRAM 11959 OF 15495 ***
    // Wavefunction(s) for diagram number 11959
    // (none)
    // Amplitude(s) for diagram number 11959
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[508], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[23] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[100] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[101] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11960 OF 15495 ***
    // Wavefunction(s) for diagram number 11960
    // (none)
    // Amplitude(s) for diagram number 11960
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[2], w_fp[500], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[484] += amp_sv[0];
    jamp_sv[485] -= amp_sv[0];
    jamp_sv[498] -= amp_sv[0];
    jamp_sv[499] += amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 11961 OF 15495 ***
    // Wavefunction(s) for diagram number 11961
    // (none)
    // Amplitude(s) for diagram number 11961
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[508], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[96] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 11962 OF 15495 ***
    // Wavefunction(s) for diagram number 11962
    // (none)
    // Amplitude(s) for diagram number 11962
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[554], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];
    jamp_sv[567] += amp_sv[0];
    jamp_sv[569] -= amp_sv[0];

    // *** DIAGRAM 11963 OF 15495 ***
    // Wavefunction(s) for diagram number 11963
    // (none)
    // Amplitude(s) for diagram number 11963
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[599], w_fp[501], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[647] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11964 OF 15495 ***
    // Wavefunction(s) for diagram number 11964
    // (none)
    // Amplitude(s) for diagram number 11964
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[597], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];

    // *** DIAGRAM 11965 OF 15495 ***
    // Wavefunction(s) for diagram number 11965
    // (none)
    // Amplitude(s) for diagram number 11965
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[696], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[33] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[57] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11966 OF 15495 ***
    // Wavefunction(s) for diagram number 11966
    // (none)
    // Amplitude(s) for diagram number 11966
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[2], w_fp[597], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[489] -= amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[495] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];

    // *** DIAGRAM 11967 OF 15495 ***
    // Wavefunction(s) for diagram number 11967
    // (none)
    // Amplitude(s) for diagram number 11967
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[696], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[31] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[55] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];

    // *** DIAGRAM 11968 OF 15495 ***
    // Wavefunction(s) for diagram number 11968
    // (none)
    // Amplitude(s) for diagram number 11968
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[599], w_fp[2], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[525] -= amp_sv[0];
    jamp_sv[527] += amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[585] += amp_sv[0];
    jamp_sv[587] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];

    // *** DIAGRAM 11969 OF 15495 ***
    // Wavefunction(s) for diagram number 11969
    // (none)
    // Amplitude(s) for diagram number 11969
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[507], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];

    // *** DIAGRAM 11970 OF 15495 ***
    // Wavefunction(s) for diagram number 11970
    // (none)
    // Amplitude(s) for diagram number 11970
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[603], w_fp[501], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];

    // *** DIAGRAM 11971 OF 15495 ***
    // Wavefunction(s) for diagram number 11971
    // (none)
    // Amplitude(s) for diagram number 11971
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[507], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[624] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[629] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11972 OF 15495 ***
    // Wavefunction(s) for diagram number 11972
    // (none)
    // Amplitude(s) for diagram number 11972
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[601], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];

    // *** DIAGRAM 11973 OF 15495 ***
    // Wavefunction(s) for diagram number 11973
    // (none)
    // Amplitude(s) for diagram number 11973
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[317], w_fp[148], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];

    // *** DIAGRAM 11974 OF 15495 ***
    // Wavefunction(s) for diagram number 11974
    // (none)
    // Amplitude(s) for diagram number 11974
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[698], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];

    // *** DIAGRAM 11975 OF 15495 ***
    // Wavefunction(s) for diagram number 11975
    // (none)
    // Amplitude(s) for diagram number 11975
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[317], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[301] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[421] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[565] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11976 OF 15495 ***
    // Wavefunction(s) for diagram number 11976
    // (none)
    // Amplitude(s) for diagram number 11976
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[2], w_fp[601], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];

    // *** DIAGRAM 11977 OF 15495 ***
    // Wavefunction(s) for diagram number 11977
    // (none)
    // Amplitude(s) for diagram number 11977
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[148], w_fp[503], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];

    // *** DIAGRAM 11978 OF 15495 ***
    // Wavefunction(s) for diagram number 11978
    // (none)
    // Amplitude(s) for diagram number 11978
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[698], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[291] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[412] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11979 OF 15495 ***
    // Wavefunction(s) for diagram number 11979
    // (none)
    // Amplitude(s) for diagram number 11979
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[2], w_fp[503], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];

    // *** DIAGRAM 11980 OF 15495 ***
    // Wavefunction(s) for diagram number 11980
    // (none)
    // Amplitude(s) for diagram number 11980
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[603], w_fp[2], w_fp[334], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[522] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[523] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[582] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11981 OF 15495 ***
    // Wavefunction(s) for diagram number 11981
    // (none)
    // Amplitude(s) for diagram number 11981
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[148], w_fp[274], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[303] -= amp_sv[0];
    jamp_sv[305] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[412] -= amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[148], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += amp_sv[0];
    jamp_sv[291] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[299] += amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[309] -= amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[148], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[300] -= amp_sv[0];
    jamp_sv[301] += amp_sv[0];
    jamp_sv[303] += amp_sv[0];
    jamp_sv[305] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[412] += amp_sv[0];
    jamp_sv[418] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];

    // *** DIAGRAM 11982 OF 15495 ***
    // Wavefunction(s) for diagram number 11982
    // (none)
    // Amplitude(s) for diagram number 11982
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[2], w_fp[274], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] -= amp_sv[0];
    jamp_sv[23] += amp_sv[0];
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[100] += amp_sv[0];
    jamp_sv[101] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[522] -= amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[582] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[2], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[76] += amp_sv[0];
    jamp_sv[77] -= amp_sv[0];
    jamp_sv[90] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];
    jamp_sv[642] -= amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[221], w_fp[2], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[22] += amp_sv[0];
    jamp_sv[23] -= amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[101] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];
    jamp_sv[484] -= amp_sv[0];
    jamp_sv[485] += amp_sv[0];
    jamp_sv[498] += amp_sv[0];
    jamp_sv[499] -= amp_sv[0];
    jamp_sv[522] += amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[582] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[628] += amp_sv[0];
    jamp_sv[629] -= amp_sv[0];

    // *** DIAGRAM 11983 OF 15495 ***
    // Wavefunction(s) for diagram number 11983
    // (none)
    // Amplitude(s) for diagram number 11983
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[519], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    jamp_sv[624] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[647] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[454], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[629] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[633] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    jamp_sv[624] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[642] += amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[647] += amp_sv[0];

    // *** DIAGRAM 11984 OF 15495 ***
    // Wavefunction(s) for diagram number 11984
    // (none)
    // Amplitude(s) for diagram number 11984
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[2], w_fp[519], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[33] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[57] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[565] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[2], w_fp[454], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[87] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[301] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[2], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[33] -= amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[57] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[291] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[489] += amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[495] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    jamp_sv[565] += amp_sv[0];

    // *** DIAGRAM 11985 OF 15495 ***
    // Wavefunction(s) for diagram number 11985
    // (none)
    // Amplitude(s) for diagram number 11985
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[554], w_fp[118], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[325] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[327] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[329] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[541] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[543] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[545] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11986 OF 15495 ***
    // Wavefunction(s) for diagram number 11986
    // (none)
    // Amplitude(s) for diagram number 11986
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[118], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[313] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[322] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[325] -= amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[329] += amp_sv[0];
    jamp_sv[332] -= amp_sv[0];
    jamp_sv[529] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[538] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[541] += amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[545] -= amp_sv[0];
    jamp_sv[548] += amp_sv[0];

    // *** DIAGRAM 11987 OF 15495 ***
    // Wavefunction(s) for diagram number 11987
    // (none)
    // Amplitude(s) for diagram number 11987
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[508], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[21] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[99] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11988 OF 15495 ***
    // Wavefunction(s) for diagram number 11988
    // (none)
    // Amplitude(s) for diagram number 11988
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[2], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[21] += amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

    // *** DIAGRAM 11989 OF 15495 ***
    // Wavefunction(s) for diagram number 11989
    // (none)
    // Amplitude(s) for diagram number 11989
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[508], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[19] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[21] -= amp_sv[0];
    jamp_sv[22] += amp_sv[0];
    jamp_sv[97] += amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[99] += amp_sv[0];
    jamp_sv[100] -= amp_sv[0];
    jamp_sv[217] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[625] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[628] += amp_sv[0];

    // *** DIAGRAM 11990 OF 15495 ***
    // Wavefunction(s) for diagram number 11990
    // (none)
    // Amplitude(s) for diagram number 11990
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[554], w_fp[2], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] -= amp_sv[0];
    jamp_sv[325] += amp_sv[0];
    jamp_sv[327] += amp_sv[0];
    jamp_sv[329] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[423] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[447] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[541] -= amp_sv[0];
    jamp_sv[543] -= amp_sv[0];
    jamp_sv[545] += amp_sv[0];

    // *** DIAGRAM 11991 OF 15495 ***
    // Wavefunction(s) for diagram number 11991
    // (none)
    // Amplitude(s) for diagram number 11991
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[607], w_fp[501], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[233] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[639] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[641] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11992 OF 15495 ***
    // Wavefunction(s) for diagram number 11992
    // (none)
    // Amplitude(s) for diagram number 11992
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[605], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[644] += amp_sv[0];

    // *** DIAGRAM 11993 OF 15495 ***
    // Wavefunction(s) for diagram number 11993
    // (none)
    // Amplitude(s) for diagram number 11993
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[246], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[39] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[81] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[315] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[531] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 11994 OF 15495 ***
    // Wavefunction(s) for diagram number 11994
    // (none)
    // Amplitude(s) for diagram number 11994
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[344], w_fp[2], w_fp[605], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[39] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[81] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[315] -= amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[531] += amp_sv[0];

    // *** DIAGRAM 11995 OF 15495 ***
    // Wavefunction(s) for diagram number 11995
    // (none)
    // Amplitude(s) for diagram number 11995
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[246], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[37] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[39] -= amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[79] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[81] += amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[313] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[315] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[529] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[531] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];

    // *** DIAGRAM 11996 OF 15495 ***
    // Wavefunction(s) for diagram number 11996
    // (none)
    // Amplitude(s) for diagram number 11996
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[607], w_fp[2], w_fp[333], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[233] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];

    // *** DIAGRAM 11997 OF 15495 ***
    // Wavefunction(s) for diagram number 11997
    // (none)
    // Amplitude(s) for diagram number 11997
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[507], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];

    // *** DIAGRAM 11998 OF 15495 ***
    // Wavefunction(s) for diagram number 11998
    // (none)
    // Amplitude(s) for diagram number 11998
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[610], w_fp[501], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];

    // *** DIAGRAM 11999 OF 15495 ***
    // Wavefunction(s) for diagram number 11999
    // (none)
    // Amplitude(s) for diagram number 11999
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[507], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[625] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[628] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12000 OF 15495 ***
    // Wavefunction(s) for diagram number 12000
    // (none)
    // Amplitude(s) for diagram number 12000
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[501], w_fp[608], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[227] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[625] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[628] -= amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[639] += amp_sv[0];
    jamp_sv[641] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];

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
