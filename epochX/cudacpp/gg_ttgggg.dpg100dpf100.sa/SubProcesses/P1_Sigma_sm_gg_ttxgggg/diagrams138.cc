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
  diagramgroup138( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 1 );
    retrieveWf( wfs, w_cx, nevt, 2 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 51 );
    retrieveWf( wfs, w_cx, nevt, 52 );
    retrieveWf( wfs, w_cx, nevt, 53 );
    retrieveWf( wfs, w_cx, nevt, 55 );
    retrieveWf( wfs, w_cx, nevt, 62 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 73 );
    retrieveWf( wfs, w_cx, nevt, 82 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 90 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 108 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 122 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 126 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 148 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 162 );
    retrieveWf( wfs, w_cx, nevt, 166 );
    retrieveWf( wfs, w_cx, nevt, 173 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 178 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 189 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 241 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 247 );
    retrieveWf( wfs, w_cx, nevt, 250 );
    retrieveWf( wfs, w_cx, nevt, 251 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 254 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 261 );
    retrieveWf( wfs, w_cx, nevt, 278 );
    retrieveWf( wfs, w_cx, nevt, 280 );
    retrieveWf( wfs, w_cx, nevt, 352 );
    retrieveWf( wfs, w_cx, nevt, 360 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 366 );
    retrieveWf( wfs, w_cx, nevt, 431 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 454 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 514 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 531 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 551 );
    retrieveWf( wfs, w_cx, nevt, 553 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 575 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 597 );
    retrieveWf( wfs, w_cx, nevt, 601 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 605 );
    retrieveWf( wfs, w_cx, nevt, 608 );
    retrieveWf( wfs, w_cx, nevt, 625 );
    retrieveWf( wfs, w_cx, nevt, 630 );
    retrieveWf( wfs, w_cx, nevt, 631 );
    retrieveWf( wfs, w_cx, nevt, 632 );
    retrieveWf( wfs, w_cx, nevt, 633 );
    retrieveWf( wfs, w_cx, nevt, 634 );
    retrieveWf( wfs, w_cx, nevt, 635 );
    retrieveWf( wfs, w_cx, nevt, 663 );
    retrieveWf( wfs, w_cx, nevt, 673 );
    retrieveWf( wfs, w_cx, nevt, 675 );
    retrieveWf( wfs, w_cx, nevt, 696 );
    retrieveWf( wfs, w_cx, nevt, 698 );
    retrieveWf( wfs, w_cx, nevt, 724 );
    retrieveWf( wfs, w_cx, nevt, 725 );
    retrieveWf( wfs, w_cx, nevt, 734 );
    retrieveWf( wfs, w_cx, nevt, 735 );
#endif
#endif

    // *** DIAGRAM 13701 OF 15495 ***
    // Wavefunction(s) for diagram number 13701
    // (none)
    // Amplitude(s) for diagram number 13701
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[595], w_fp[2], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[464] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];

    // *** DIAGRAM 13702 OF 15495 ***
    // Wavefunction(s) for diagram number 13702
    // (none)
    // Amplitude(s) for diagram number 13702
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[173], w_fp[122], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13703 OF 15495 ***
    // Wavefunction(s) for diagram number 13703
    // (none)
    // Amplitude(s) for diagram number 13703
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[595], w_fp[122], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[464] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13704 OF 15495 ***
    // Wavefunction(s) for diagram number 13704
    // (none)
    // Amplitude(s) for diagram number 13704
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[173], w_fp[2], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13705 OF 15495 ***
    // Wavefunction(s) for diagram number 13705
    // (none)
    // Amplitude(s) for diagram number 13705
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[124], w_fp[724], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13706 OF 15495 ***
    // Wavefunction(s) for diagram number 13706
    // (none)
    // Amplitude(s) for diagram number 13706
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[261], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[352], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[278], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[284] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[464] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13707 OF 15495 ***
    // Wavefunction(s) for diagram number 13707
    // (none)
    // Amplitude(s) for diagram number 13707
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[663], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 13708 OF 15495 ***
    // Wavefunction(s) for diagram number 13708
    // (none)
    // Amplitude(s) for diagram number 13708
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[453], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13709 OF 15495 ***
    // Wavefunction(s) for diagram number 13709
    // (none)
    // Amplitude(s) for diagram number 13709
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[522], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 13710 OF 15495 ***
    // Wavefunction(s) for diagram number 13710
    // (none)
    // Amplitude(s) for diagram number 13710
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[4], w_fp[673], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[268] -= amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 13711 OF 15495 ***
    // Wavefunction(s) for diagram number 13711
    // (none)
    // Amplitude(s) for diagram number 13711
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[150], w_fp[453], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[148] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[458] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 13712 OF 15495 ***
    // Wavefunction(s) for diagram number 13712
    // (none)
    // Amplitude(s) for diagram number 13712
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[673], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13713 OF 15495 ***
    // Wavefunction(s) for diagram number 13713
    // (none)
    // Amplitude(s) for diagram number 13713
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[663], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 13714 OF 15495 ***
    // Wavefunction(s) for diagram number 13714
    // (none)
    // Amplitude(s) for diagram number 13714
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[189], w_fp[122], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13715 OF 15495 ***
    // Wavefunction(s) for diagram number 13715
    // (none)
    // Amplitude(s) for diagram number 13715
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[675], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13716 OF 15495 ***
    // Wavefunction(s) for diagram number 13716
    // (none)
    // Amplitude(s) for diagram number 13716
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[189], w_fp[2], w_fp[124], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13717 OF 15495 ***
    // Wavefunction(s) for diagram number 13717
    // (none)
    // Amplitude(s) for diagram number 13717
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[625], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13718 OF 15495 ***
    // Wavefunction(s) for diagram number 13718
    // (none)
    // Amplitude(s) for diagram number 13718
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[360], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 13719 OF 15495 ***
    // Wavefunction(s) for diagram number 13719
    // (none)
    // Amplitude(s) for diagram number 13719
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[4], w_fp[55], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 13720 OF 15495 ***
    // Wavefunction(s) for diagram number 13720
    // (none)
    // Amplitude(s) for diagram number 13720
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[360], w_fp[4], w_fp[725], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 13721 OF 15495 ***
    // Wavefunction(s) for diagram number 13721
    // (none)
    // Amplitude(s) for diagram number 13721
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[55], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13722 OF 15495 ***
    // Wavefunction(s) for diagram number 13722
    // (none)
    // Amplitude(s) for diagram number 13722
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[2], w_fp[360], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 13723 OF 15495 ***
    // Wavefunction(s) for diagram number 13723
    // (none)
    // Amplitude(s) for diagram number 13723
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[150], w_fp[124], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[150], w_fp[124], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[150], w_fp[124], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13724 OF 15495 ***
    // Wavefunction(s) for diagram number 13724
    // (none)
    // Amplitude(s) for diagram number 13724
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[124], w_fp[725], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13725 OF 15495 ***
    // Wavefunction(s) for diagram number 13725
    // (none)
    // Amplitude(s) for diagram number 13725
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[150], w_fp[625], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];

    // *** DIAGRAM 13726 OF 15495 ***
    // Wavefunction(s) for diagram number 13726
    // (none)
    // Amplitude(s) for diagram number 13726
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[675], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13727 OF 15495 ***
    // Wavefunction(s) for diagram number 13727
    // (none)
    // Amplitude(s) for diagram number 13727
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[122], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13728 OF 15495 ***
    // Wavefunction(s) for diagram number 13728
    // (none)
    // Amplitude(s) for diagram number 13728
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[90], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[244] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[126], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[62], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[244] -= amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[268] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 13729 OF 15495 ***
    // Wavefunction(s) for diagram number 13729
    // (none)
    // Amplitude(s) for diagram number 13729
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[90], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[126], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[62], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[20] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13730 OF 15495 ***
    // Wavefunction(s) for diagram number 13730
    // (none)
    // Amplitude(s) for diagram number 13730
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[475], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[553], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[592], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[260] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[458] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13731 OF 15495 ***
    // Wavefunction(s) for diagram number 13731
    // (none)
    // Amplitude(s) for diagram number 13731
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[475], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[553], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[458] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[592], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 13732 OF 15495 ***
    // Wavefunction(s) for diagram number 13732
    // (none)
    // Amplitude(s) for diagram number 13732
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[261], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[352], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[278], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[4] -= amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[28] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[464] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13733 OF 15495 ***
    // Wavefunction(s) for diagram number 13733
    // (none)
    // Amplitude(s) for diagram number 13733
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[108], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[431], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[280], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[308] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[686] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13734 OF 15495 ***
    // Wavefunction(s) for diagram number 13734
    // (none)
    // Amplitude(s) for diagram number 13734
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[51], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[166], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[241], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];

    // *** DIAGRAM 13735 OF 15495 ***
    // Wavefunction(s) for diagram number 13735
    // (none)
    // Amplitude(s) for diagram number 13735
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[591], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[539], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[2], w_fp[531], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[110] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13736 OF 15495 ***
    // Wavefunction(s) for diagram number 13736
    // (none)
    // Amplitude(s) for diagram number 13736
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[591], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[539], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[68] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[164] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[531], w_fp[1], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13737 OF 15495 ***
    // Wavefunction(s) for diagram number 13737
    // (none)
    // Amplitude(s) for diagram number 13737
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[51], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[166], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[255], w_fp[241], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];

    // *** DIAGRAM 13738 OF 15495 ***
    // Wavefunction(s) for diagram number 13738
    // (none)
    // Amplitude(s) for diagram number 13738
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[108], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[431], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[280], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];

    // *** DIAGRAM 13739 OF 15495 ***
    // Wavefunction(s) for diagram number 13739
    // (none)
    // Amplitude(s) for diagram number 13739
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[52], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] -= amp_sv[0];
    jamp_sv[4] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[20] -= amp_sv[0];
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[308] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[686] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[53], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[110] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[306] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[684] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[162], COUPs[1], 1.0, &amp_fp[0] );
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
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[635], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[4] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[634], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[633], COUPs[1], 1.0, &amp_fp[0] );
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
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[632], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[1] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[18] -= amp_sv[0];
    jamp_sv[20] += amp_sv[0];
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[308] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[686] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[631], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[34] -= amp_sv[0];
    jamp_sv[58] += amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[110] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[306] += amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[684] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[630], COUPs[1], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 13740 OF 15495 ***
    // Wavefunction(s) for diagram number 13740
    // (none)
    // Amplitude(s) for diagram number 13740
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[254], w_fp[6], w_fp[735], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13741 OF 15495 ***
    // Wavefunction(s) for diagram number 13741
    // (none)
    // Amplitude(s) for diagram number 13741
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[2], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[416] += amp_sv[0];

    // *** DIAGRAM 13742 OF 15495 ***
    // Wavefunction(s) for diagram number 13742
    // (none)
    // Amplitude(s) for diagram number 13742
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[178], w_fp[148], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13743 OF 15495 ***
    // Wavefunction(s) for diagram number 13743
    // (none)
    // Amplitude(s) for diagram number 13743
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[148], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13744 OF 15495 ***
    // Wavefunction(s) for diagram number 13744
    // (none)
    // Amplitude(s) for diagram number 13744
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[178], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[302] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13745 OF 15495 ***
    // Wavefunction(s) for diagram number 13745
    // (none)
    // Amplitude(s) for diagram number 13745
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[68], w_fp[735], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13746 OF 15495 ***
    // Wavefunction(s) for diagram number 13746
    // (none)
    // Amplitude(s) for diagram number 13746
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[247], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[251], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[250], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[512] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[566] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13747 OF 15495 ***
    // Wavefunction(s) for diagram number 13747
    // (none)
    // Amplitude(s) for diagram number 13747
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[696], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];

    // *** DIAGRAM 13748 OF 15495 ***
    // Wavefunction(s) for diagram number 13748
    // (none)
    // Amplitude(s) for diagram number 13748
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[597], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13749 OF 15495 ***
    // Wavefunction(s) for diagram number 13749
    // (none)
    // Amplitude(s) for diagram number 13749
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[514], w_fp[1], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13750 OF 15495 ***
    // Wavefunction(s) for diagram number 13750
    // (none)
    // Amplitude(s) for diagram number 13750
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[6], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[504] += amp_sv[0];
    jamp_sv[506] -= amp_sv[0];
    jamp_sv[512] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 13751 OF 15495 ***
    // Wavefunction(s) for diagram number 13751
    // (none)
    // Amplitude(s) for diagram number 13751
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[144], w_fp[597], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[290] += amp_sv[0];
    jamp_sv[410] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13752 OF 15495 ***
    // Wavefunction(s) for diagram number 13752
    // (none)
    // Amplitude(s) for diagram number 13752
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13753 OF 15495 ***
    // Wavefunction(s) for diagram number 13753
    // (none)
    // Amplitude(s) for diagram number 13753
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[696], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += amp_sv[0];
    jamp_sv[54] -= amp_sv[0];
    jamp_sv[288] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];

    // *** DIAGRAM 13754 OF 15495 ***
    // Wavefunction(s) for diagram number 13754
    // (none)
    // Amplitude(s) for diagram number 13754
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[187], w_fp[148], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13755 OF 15495 ***
    // Wavefunction(s) for diagram number 13755
    // (none)
    // Amplitude(s) for diagram number 13755
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[698], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13756 OF 15495 ***
    // Wavefunction(s) for diagram number 13756
    // (none)
    // Amplitude(s) for diagram number 13756
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[187], w_fp[2], w_fp[68], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13757 OF 15495 ***
    // Wavefunction(s) for diagram number 13757
    // (none)
    // Amplitude(s) for diagram number 13757
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[601], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13758 OF 15495 ***
    // Wavefunction(s) for diagram number 13758
    // (none)
    // Amplitude(s) for diagram number 13758
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[254], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];

    // *** DIAGRAM 13759 OF 15495 ***
    // Wavefunction(s) for diagram number 13759
    // (none)
    // Amplitude(s) for diagram number 13759
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[6], w_fp[182], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13760 OF 15495 ***
    // Wavefunction(s) for diagram number 13760
    // (none)
    // Amplitude(s) for diagram number 13760
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[254], w_fp[6], w_fp[734], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13761 OF 15495 ***
    // Wavefunction(s) for diagram number 13761
    // (none)
    // Amplitude(s) for diagram number 13761
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[182], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13762 OF 15495 ***
    // Wavefunction(s) for diagram number 13762
    // (none)
    // Amplitude(s) for diagram number 13762
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[2], w_fp[254], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];

    // *** DIAGRAM 13763 OF 15495 ***
    // Wavefunction(s) for diagram number 13763
    // (none)
    // Amplitude(s) for diagram number 13763
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[144], w_fp[68], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[144], w_fp[68], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[144], w_fp[68], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13764 OF 15495 ***
    // Wavefunction(s) for diagram number 13764
    // (none)
    // Amplitude(s) for diagram number 13764
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[68], w_fp[734], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[14] += amp_sv[0];
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[302] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13765 OF 15495 ***
    // Wavefunction(s) for diagram number 13765
    // (none)
    // Amplitude(s) for diagram number 13765
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[144], w_fp[601], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];

    // *** DIAGRAM 13766 OF 15495 ***
    // Wavefunction(s) for diagram number 13766
    // (none)
    // Amplitude(s) for diagram number 13766
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[698], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13767 OF 15495 ***
    // Wavefunction(s) for diagram number 13767
    // (none)
    // Amplitude(s) for diagram number 13767
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[148], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13768 OF 15495 ***
    // Wavefunction(s) for diagram number 13768
    // (none)
    // Amplitude(s) for diagram number 13768
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[136], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[480] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[137], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[494] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[551], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[288] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[480] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[504] -= amp_sv[0];
    jamp_sv[506] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];

    // *** DIAGRAM 13769 OF 15495 ***
    // Wavefunction(s) for diagram number 13769
    // (none)
    // Amplitude(s) for diagram number 13769
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[136], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[137], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[30] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[54] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[551], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13770 OF 15495 ***
    // Wavefunction(s) for diagram number 13770
    // (none)
    // Amplitude(s) for diagram number 13770
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[519], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[454], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[86] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[300] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[575], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[488] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[494] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[564] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13771 OF 15495 ***
    // Wavefunction(s) for diagram number 13771
    // (none)
    // Amplitude(s) for diagram number 13771
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[519], w_fp[1], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[540] -= amp_sv[0];
    jamp_sv[564] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[454], w_fp[1], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[300] += amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[575], w_fp[1], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[290] -= amp_sv[0];
    jamp_sv[410] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[494] += amp_sv[0];
    jamp_sv[540] += amp_sv[0];
    jamp_sv[564] -= amp_sv[0];

    // *** DIAGRAM 13772 OF 15495 ***
    // Wavefunction(s) for diagram number 13772
    // (none)
    // Amplitude(s) for diagram number 13772
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[247], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] -= amp_sv[0];
    jamp_sv[6] += amp_sv[0];
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[86] += amp_sv[0];
    jamp_sv[152] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[566] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[251], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += amp_sv[0];
    jamp_sv[14] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[32] += amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[56] -= amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[302] -= amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[250], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[0] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[54] += amp_sv[0];
    jamp_sv[72] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[86] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[416] -= amp_sv[0];
    jamp_sv[512] += amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[566] += amp_sv[0];

    // *** DIAGRAM 13773 OF 15495 ***
    // Wavefunction(s) for diagram number 13773
    // (none)
    // Amplitude(s) for diagram number 13773
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[362], w_fp[5], w_fp[735], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13774 OF 15495 ***
    // Wavefunction(s) for diagram number 13774
    // (none)
    // Amplitude(s) for diagram number 13774
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[2], w_fp[362], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[536] += amp_sv[0];

    // *** DIAGRAM 13775 OF 15495 ***
    // Wavefunction(s) for diagram number 13775
    // (none)
    // Amplitude(s) for diagram number 13775
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[178], w_fp[118], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13776 OF 15495 ***
    // Wavefunction(s) for diagram number 13776
    // (none)
    // Amplitude(s) for diagram number 13776
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[438], w_fp[118], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13777 OF 15495 ***
    // Wavefunction(s) for diagram number 13777
    // (none)
    // Amplitude(s) for diagram number 13777
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[178], w_fp[2], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[326] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 13778 OF 15495 ***
    // Wavefunction(s) for diagram number 13778
    // (none)
    // Amplitude(s) for diagram number 13778
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[87], w_fp[735], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13779 OF 15495 ***
    // Wavefunction(s) for diagram number 13779
    // (none)
    // Amplitude(s) for diagram number 13779
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[82], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[73], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[366], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13780 OF 15495 ***
    // Wavefunction(s) for diagram number 13780
    // (none)
    // Amplitude(s) for diagram number 13780
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[246], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];

    // *** DIAGRAM 13781 OF 15495 ***
    // Wavefunction(s) for diagram number 13781
    // (none)
    // Amplitude(s) for diagram number 13781
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[605], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13782 OF 15495 ***
    // Wavefunction(s) for diagram number 13782
    // (none)
    // Amplitude(s) for diagram number 13782
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[1], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 13783 OF 15495 ***
    // Wavefunction(s) for diagram number 13783
    // (none)
    // Amplitude(s) for diagram number 13783
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[5], w_fp[604], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 13784 OF 15495 ***
    // Wavefunction(s) for diagram number 13784
    // (none)
    // Amplitude(s) for diagram number 13784
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[144], w_fp[605], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];

    // *** DIAGRAM 13785 OF 15495 ***
    // Wavefunction(s) for diagram number 13785
    // (none)
    // Amplitude(s) for diagram number 13785
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[2], w_fp[604], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13786 OF 15495 ***
    // Wavefunction(s) for diagram number 13786
    // (none)
    // Amplitude(s) for diagram number 13786
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[246], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] += amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];

    // *** DIAGRAM 13787 OF 15495 ***
    // Wavefunction(s) for diagram number 13787
    // (none)
    // Amplitude(s) for diagram number 13787
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[187], w_fp[118], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13788 OF 15495 ***
    // Wavefunction(s) for diagram number 13788
    // (none)
    // Amplitude(s) for diagram number 13788
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[92], w_fp[5], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13789 OF 15495 ***
    // Wavefunction(s) for diagram number 13789
    // (none)
    // Amplitude(s) for diagram number 13789
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[187], w_fp[2], w_fp[87], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[324] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13790 OF 15495 ***
    // Wavefunction(s) for diagram number 13790
    // (none)
    // Amplitude(s) for diagram number 13790
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[252], w_fp[2], w_fp[608], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[56] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[62] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[324] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[540] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13791 OF 15495 ***
    // Wavefunction(s) for diagram number 13791
    // (none)
    // Amplitude(s) for diagram number 13791
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[362], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[362], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[362], w_fp[144], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[156] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[534] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];

    // *** DIAGRAM 13792 OF 15495 ***
    // Wavefunction(s) for diagram number 13792
    // (none)
    // Amplitude(s) for diagram number 13792
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[5], w_fp[64], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[318] += amp_sv[0];
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[534] -= amp_sv[0];

    // *** DIAGRAM 13793 OF 15495 ***
    // Wavefunction(s) for diagram number 13793
    // (none)
    // Amplitude(s) for diagram number 13793
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[362], w_fp[5], w_fp[734], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];

    // *** DIAGRAM 13794 OF 15495 ***
    // Wavefunction(s) for diagram number 13794
    // (none)
    // Amplitude(s) for diagram number 13794
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[2], w_fp[64], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[36] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13795 OF 15495 ***
    // Wavefunction(s) for diagram number 13795
    // (none)
    // Amplitude(s) for diagram number 13795
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[2], w_fp[362], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += amp_sv[0];
    jamp_sv[198] -= amp_sv[0];
    jamp_sv[318] -= amp_sv[0];
    jamp_sv[534] += amp_sv[0];

    // *** DIAGRAM 13796 OF 15495 ***
    // Wavefunction(s) for diagram number 13796
    // (none)
    // Amplitude(s) for diagram number 13796
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[144], w_fp[87], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[144], w_fp[87], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[1], w_fp[144], w_fp[87], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13797 OF 15495 ***
    // Wavefunction(s) for diagram number 13797
    // (none)
    // Amplitude(s) for diagram number 13797
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[87], w_fp[734], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[2] += amp_sv[0];
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[12] -= amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 13798 OF 15495 ***
    // Wavefunction(s) for diagram number 13798
    // (none)
    // Amplitude(s) for diagram number 13798
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[1], w_fp[144], w_fp[608], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[56] += amp_sv[0];
    jamp_sv[62] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[324] += amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[540] -= amp_sv[0];

    // *** DIAGRAM 13799 OF 15495 ***
    // Wavefunction(s) for diagram number 13799
    // (none)
    // Amplitude(s) for diagram number 13799
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[180], w_fp[92], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 13800 OF 15495 ***
    // Wavefunction(s) for diagram number 13800
    // (none)
    // Amplitude(s) for diagram number 13800
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[557], w_fp[118], w_fp[1], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[318] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[534] -= cxtype( 0, 1 ) * amp_sv[0];

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
