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
  diagramgroup109( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 11 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 13 );
    retrieveWf( wfs, w_cx, nevt, 16 );
    retrieveWf( wfs, w_cx, nevt, 27 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 47 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 92 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 110 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 242 );
    retrieveWf( wfs, w_cx, nevt, 245 );
    retrieveWf( wfs, w_cx, nevt, 246 );
    retrieveWf( wfs, w_cx, nevt, 267 );
    retrieveWf( wfs, w_cx, nevt, 268 );
    retrieveWf( wfs, w_cx, nevt, 274 );
    retrieveWf( wfs, w_cx, nevt, 285 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 302 );
    retrieveWf( wfs, w_cx, nevt, 304 );
    retrieveWf( wfs, w_cx, nevt, 309 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 493 );
    retrieveWf( wfs, w_cx, nevt, 494 );
    retrieveWf( wfs, w_cx, nevt, 497 );
    retrieveWf( wfs, w_cx, nevt, 498 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 548 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 593 );
    retrieveWf( wfs, w_cx, nevt, 606 );
    retrieveWf( wfs, w_cx, nevt, 607 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 610 );
    retrieveWf( wfs, w_cx, nevt, 613 );
    retrieveWf( wfs, w_cx, nevt, 614 );
    retrieveWf( wfs, w_cx, nevt, 616 );
    retrieveWf( wfs, w_cx, nevt, 617 );
    retrieveWf( wfs, w_cx, nevt, 628 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 670 );
    retrieveWf( wfs, w_cx, nevt, 697 );
    retrieveWf( wfs, w_cx, nevt, 703 );
    retrieveWf( wfs, w_cx, nevt, 706 );
    retrieveWf( wfs, w_cx, nevt, 707 );
    retrieveWf( wfs, w_cx, nevt, 708 );
    retrieveWf( wfs, w_cx, nevt, 709 );
    retrieveWf( wfs, w_cx, nevt, 710 );
    retrieveWf( wfs, w_cx, nevt, 711 );
    retrieveWf( wfs, w_cx, nevt, 712 );
    retrieveWf( wfs, w_cx, nevt, 713 );
    retrieveWf( wfs, w_cx, nevt, 714 );
    retrieveWf( wfs, w_cx, nevt, 743 );
    retrieveWf( wfs, w_cx, nevt, 748 );
#endif
#endif

    // *** DIAGRAM 10801 OF 15495 ***
    // Wavefunction(s) for diagram number 10801
    // (none)
    // Amplitude(s) for diagram number 10801
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[309], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];

    // *** DIAGRAM 10802 OF 15495 ***
    // Wavefunction(s) for diagram number 10802
    // (none)
    // Amplitude(s) for diagram number 10802
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[2], w_fp[285], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[78] -= amp_sv[0];
    jamp_sv[80] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[2], w_fp[268], COUPs[1], 1.0, &amp_fp[0] );
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
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[2], w_fp[267], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];

    // *** DIAGRAM 10803 OF 15495 ***
    // Wavefunction(s) for diagram number 10803
    // (none)
    // Amplitude(s) for diagram number 10803
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[548], w_fp[118], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[333] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[335] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[549] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[551] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10804 OF 15495 ***
    // Wavefunction(s) for diagram number 10804
    // (none)
    // Amplitude(s) for diagram number 10804
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[118], w_fp[707], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];

    // *** DIAGRAM 10805 OF 15495 ***
    // Wavefunction(s) for diagram number 10805
    // (none)
    // Amplitude(s) for diagram number 10805
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[500], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[11] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10806 OF 15495 ***
    // Wavefunction(s) for diagram number 10806
    // (none)
    // Amplitude(s) for diagram number 10806
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[2], w_fp[707], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[603] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 10807 OF 15495 ***
    // Wavefunction(s) for diagram number 10807
    // (none)
    // Amplitude(s) for diagram number 10807
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[500], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];

    // *** DIAGRAM 10808 OF 15495 ***
    // Wavefunction(s) for diagram number 10808
    // (none)
    // Amplitude(s) for diagram number 10808
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[548], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 10809 OF 15495 ***
    // Wavefunction(s) for diagram number 10809
    // (none)
    // Amplitude(s) for diagram number 10809
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[607], w_fp[494], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[189] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[191] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10810 OF 15495 ***
    // Wavefunction(s) for diagram number 10810
    // (none)
    // Amplitude(s) for diagram number 10810
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[606], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];

    // *** DIAGRAM 10811 OF 15495 ***
    // Wavefunction(s) for diagram number 10811
    // (none)
    // Amplitude(s) for diagram number 10811
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[246], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[41] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[82] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[83] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10812 OF 15495 ***
    // Wavefunction(s) for diagram number 10812
    // (none)
    // Amplitude(s) for diagram number 10812
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[606], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[610] += amp_sv[0];
    jamp_sv[611] -= amp_sv[0];
    jamp_sv[620] -= amp_sv[0];
    jamp_sv[621] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 10813 OF 15495 ***
    // Wavefunction(s) for diagram number 10813
    // (none)
    // Amplitude(s) for diagram number 10813
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[246], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];

    // *** DIAGRAM 10814 OF 15495 ***
    // Wavefunction(s) for diagram number 10814
    // (none)
    // Amplitude(s) for diagram number 10814
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[607], w_fp[2], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[639] -= amp_sv[0];
    jamp_sv[641] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];
    jamp_sv[681] += amp_sv[0];
    jamp_sv[683] -= amp_sv[0];

    // *** DIAGRAM 10815 OF 15495 ***
    // Wavefunction(s) for diagram number 10815
    // (none)
    // Amplitude(s) for diagram number 10815
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[498], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 10816 OF 15495 ***
    // Wavefunction(s) for diagram number 10816
    // (none)
    // Amplitude(s) for diagram number 10816
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[610], w_fp[494], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];

    // *** DIAGRAM 10817 OF 15495 ***
    // Wavefunction(s) for diagram number 10817
    // (none)
    // Amplitude(s) for diagram number 10817
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[498], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10818 OF 15495 ***
    // Wavefunction(s) for diagram number 10818
    // (none)
    // Amplitude(s) for diagram number 10818
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[609], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];

    // *** DIAGRAM 10819 OF 15495 ***
    // Wavefunction(s) for diagram number 10819
    // (none)
    // Amplitude(s) for diagram number 10819
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[118], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];

    // *** DIAGRAM 10820 OF 15495 ***
    // Wavefunction(s) for diagram number 10820
    // (none)
    // Amplitude(s) for diagram number 10820
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[92], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];

    // *** DIAGRAM 10821 OF 15495 ***
    // Wavefunction(s) for diagram number 10821
    // (none)
    // Amplitude(s) for diagram number 10821
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[2], w_fp[88], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10822 OF 15495 ***
    // Wavefunction(s) for diagram number 10822
    // (none)
    // Amplitude(s) for diagram number 10822
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[609], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 10823 OF 15495 ***
    // Wavefunction(s) for diagram number 10823
    // (none)
    // Amplitude(s) for diagram number 10823
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[118], w_fp[493], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];

    // *** DIAGRAM 10824 OF 15495 ***
    // Wavefunction(s) for diagram number 10824
    // (none)
    // Amplitude(s) for diagram number 10824
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[92], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[316] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[317] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[532] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[533] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10825 OF 15495 ***
    // Wavefunction(s) for diagram number 10825
    // (none)
    // Amplitude(s) for diagram number 10825
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[2], w_fp[493], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];

    // *** DIAGRAM 10826 OF 15495 ***
    // Wavefunction(s) for diagram number 10826
    // (none)
    // Amplitude(s) for diagram number 10826
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[610], w_fp[2], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[187] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[679] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10827 OF 15495 ***
    // Wavefunction(s) for diagram number 10827
    // (none)
    // Amplitude(s) for diagram number 10827
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[118], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[118], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[316] += amp_sv[0];
    jamp_sv[317] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[321] += amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[327] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[532] -= amp_sv[0];
    jamp_sv[533] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[537] -= amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[543] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[118], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];

    // *** DIAGRAM 10828 OF 15495 ***
    // Wavefunction(s) for diagram number 10828
    // (none)
    // Amplitude(s) for diagram number 10828
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[2], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[679] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[2], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += amp_sv[0];
    jamp_sv[99] -= amp_sv[0];
    jamp_sv[108] -= amp_sv[0];
    jamp_sv[109] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[206], w_fp[2], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[603] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    jamp_sv[679] += amp_sv[0];

    // *** DIAGRAM 10829 OF 15495 ***
    // Wavefunction(s) for diagram number 10829
    // (none)
    // Amplitude(s) for diagram number 10829
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[558], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 10830 OF 15495 ***
    // Wavefunction(s) for diagram number 10830
    // (none)
    // Amplitude(s) for diagram number 10830
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[550], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[480], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[558], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[316] -= amp_sv[0];
    jamp_sv[317] += amp_sv[0];
    jamp_sv[532] += amp_sv[0];
    jamp_sv[533] -= amp_sv[0];
    jamp_sv[610] -= amp_sv[0];
    jamp_sv[611] += amp_sv[0];
    jamp_sv[620] += amp_sv[0];
    jamp_sv[621] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 10831 OF 15495 ***
    // Wavefunction(s) for diagram number 10831
    // (none)
    // Amplitude(s) for diagram number 10831
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[548], w_fp[98], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[357] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[359] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[671] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10832 OF 15495 ***
    // Wavefunction(s) for diagram number 10832
    // (none)
    // Amplitude(s) for diagram number 10832
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], w_fp[708], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 10833 OF 15495 ***
    // Wavefunction(s) for diagram number 10833
    // (none)
    // Amplitude(s) for diagram number 10833
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[500], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[9] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10834 OF 15495 ***
    // Wavefunction(s) for diagram number 10834
    // (none)
    // Amplitude(s) for diagram number 10834
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[2], w_fp[708], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[483] -= amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];

    // *** DIAGRAM 10835 OF 15495 ***
    // Wavefunction(s) for diagram number 10835
    // (none)
    // Amplitude(s) for diagram number 10835
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[500], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];

    // *** DIAGRAM 10836 OF 15495 ***
    // Wavefunction(s) for diagram number 10836
    // (none)
    // Amplitude(s) for diagram number 10836
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[548], w_fp[2], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[549] -= amp_sv[0];
    jamp_sv[551] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 10837 OF 15495 ***
    // Wavefunction(s) for diagram number 10837
    // (none)
    // Amplitude(s) for diagram number 10837
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[614], w_fp[494], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[183] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[185] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10838 OF 15495 ***
    // Wavefunction(s) for diagram number 10838
    // (none)
    // Amplitude(s) for diagram number 10838
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[613], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];

    // *** DIAGRAM 10839 OF 15495 ***
    // Wavefunction(s) for diagram number 10839
    // (none)
    // Amplitude(s) for diagram number 10839
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[245], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[47] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[106] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[107] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10840 OF 15495 ***
    // Wavefunction(s) for diagram number 10840
    // (none)
    // Amplitude(s) for diagram number 10840
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[613], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[490] += amp_sv[0];
    jamp_sv[491] -= amp_sv[0];
    jamp_sv[500] -= amp_sv[0];
    jamp_sv[501] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 10841 OF 15495 ***
    // Wavefunction(s) for diagram number 10841
    // (none)
    // Amplitude(s) for diagram number 10841
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[245], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];

    // *** DIAGRAM 10842 OF 15495 ***
    // Wavefunction(s) for diagram number 10842
    // (none)
    // Amplitude(s) for diagram number 10842
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[614], w_fp[2], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[519] -= amp_sv[0];
    jamp_sv[521] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];
    jamp_sv[561] += amp_sv[0];
    jamp_sv[563] -= amp_sv[0];

    // *** DIAGRAM 10843 OF 15495 ***
    // Wavefunction(s) for diagram number 10843
    // (none)
    // Amplitude(s) for diagram number 10843
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[498], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];

    // *** DIAGRAM 10844 OF 15495 ***
    // Wavefunction(s) for diagram number 10844
    // (none)
    // Amplitude(s) for diagram number 10844
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[617], w_fp[494], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];

    // *** DIAGRAM 10845 OF 15495 ***
    // Wavefunction(s) for diagram number 10845
    // (none)
    // Amplitude(s) for diagram number 10845
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[498], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10846 OF 15495 ***
    // Wavefunction(s) for diagram number 10846
    // (none)
    // Amplitude(s) for diagram number 10846
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[616], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];

    // *** DIAGRAM 10847 OF 15495 ***
    // Wavefunction(s) for diagram number 10847
    // (none)
    // Amplitude(s) for diagram number 10847
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[98], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];

    // *** DIAGRAM 10848 OF 15495 ***
    // Wavefunction(s) for diagram number 10848
    // (none)
    // Amplitude(s) for diagram number 10848
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[110], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[340] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];

    // *** DIAGRAM 10849 OF 15495 ***
    // Wavefunction(s) for diagram number 10849
    // (none)
    // Amplitude(s) for diagram number 10849
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[2], w_fp[104], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10850 OF 15495 ***
    // Wavefunction(s) for diagram number 10850
    // (none)
    // Amplitude(s) for diagram number 10850
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[616], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];

    // *** DIAGRAM 10851 OF 15495 ***
    // Wavefunction(s) for diagram number 10851
    // (none)
    // Amplitude(s) for diagram number 10851
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], w_fp[748], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];

    // *** DIAGRAM 10852 OF 15495 ***
    // Wavefunction(s) for diagram number 10852
    // (none)
    // Amplitude(s) for diagram number 10852
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[110], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[341] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[653] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10853 OF 15495 ***
    // Wavefunction(s) for diagram number 10853
    // (none)
    // Amplitude(s) for diagram number 10853
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[2], w_fp[748], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];

    // *** DIAGRAM 10854 OF 15495 ***
    // Wavefunction(s) for diagram number 10854
    // (none)
    // Amplitude(s) for diagram number 10854
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[617], w_fp[2], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[181] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[559] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10855 OF 15495 ***
    // Wavefunction(s) for diagram number 10855
    // (none)
    // Amplitude(s) for diagram number 10855
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[340] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], w_fp[703], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[340] += amp_sv[0];
    jamp_sv[341] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[345] += amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[351] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[652] -= amp_sv[0];
    jamp_sv[653] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[657] -= amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[663] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[98], w_fp[670], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[357] += amp_sv[0];
    jamp_sv[359] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[671] += amp_sv[0];

    // *** DIAGRAM 10856 OF 15495 ***
    // Wavefunction(s) for diagram number 10856
    // (none)
    // Amplitude(s) for diagram number 10856
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[2], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[559] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[2], w_fp[703], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += amp_sv[0];
    jamp_sv[75] -= amp_sv[0];
    jamp_sv[84] -= amp_sv[0];
    jamp_sv[85] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[181] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[208], w_fp[2], w_fp[670], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[483] += amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    jamp_sv[559] += amp_sv[0];

    // *** DIAGRAM 10857 OF 15495 ***
    // Wavefunction(s) for diagram number 10857
    // (none)
    // Amplitude(s) for diagram number 10857
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[183] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[529], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[181] += amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];

    // *** DIAGRAM 10858 OF 15495 ***
    // Wavefunction(s) for diagram number 10858
    // (none)
    // Amplitude(s) for diagram number 10858
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[486], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[106] += amp_sv[0];
    jamp_sv[107] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[547], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[82] += amp_sv[0];
    jamp_sv[83] -= amp_sv[0];
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[529], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[340] -= amp_sv[0];
    jamp_sv[341] += amp_sv[0];
    jamp_sv[490] -= amp_sv[0];
    jamp_sv[491] += amp_sv[0];
    jamp_sv[500] += amp_sv[0];
    jamp_sv[501] -= amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[652] += amp_sv[0];
    jamp_sv[653] -= amp_sv[0];

    // *** DIAGRAM 10859 OF 15495 ***
    // Wavefunction(s) for diagram number 10859
    // (none)
    // Amplitude(s) for diagram number 10859
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[548], w_fp[128], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[591] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[593] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[713] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10860 OF 15495 ***
    // Wavefunction(s) for diagram number 10860
    // (none)
    // Amplitude(s) for diagram number 10860
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[128], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 10861 OF 15495 ***
    // Wavefunction(s) for diagram number 10861
    // (none)
    // Amplitude(s) for diagram number 10861
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[500], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10862 OF 15495 ***
    // Wavefunction(s) for diagram number 10862
    // (none)
    // Amplitude(s) for diagram number 10862
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[2], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];

    // *** DIAGRAM 10863 OF 15495 ***
    // Wavefunction(s) for diagram number 10863
    // (none)
    // Amplitude(s) for diagram number 10863
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[500], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];

    // *** DIAGRAM 10864 OF 15495 ***
    // Wavefunction(s) for diagram number 10864
    // (none)
    // Amplitude(s) for diagram number 10864
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[548], w_fp[2], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];

    // *** DIAGRAM 10865 OF 15495 ***
    // Wavefunction(s) for diagram number 10865
    // (none)
    // Amplitude(s) for diagram number 10865
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10866 OF 15495 ***
    // Wavefunction(s) for diagram number 10866
    // (none)
    // Amplitude(s) for diagram number 10866
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];

    // *** DIAGRAM 10867 OF 15495 ***
    // Wavefunction(s) for diagram number 10867
    // (none)
    // Amplitude(s) for diagram number 10867
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[497], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[93] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[116] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[117] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10868 OF 15495 ***
    // Wavefunction(s) for diagram number 10868
    // (none)
    // Amplitude(s) for diagram number 10868
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[546], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] -= amp_sv[0];
    jamp_sv[93] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[256] += amp_sv[0];
    jamp_sv[257] -= amp_sv[0];
    jamp_sv[262] -= amp_sv[0];
    jamp_sv[263] += amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 10869 OF 15495 ***
    // Wavefunction(s) for diagram number 10869
    // (none)
    // Amplitude(s) for diagram number 10869
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[497], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 10870 OF 15495 ***
    // Wavefunction(s) for diagram number 10870
    // (none)
    // Amplitude(s) for diagram number 10870
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[442], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[273] -= amp_sv[0];
    jamp_sv[275] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[297] += amp_sv[0];
    jamp_sv[299] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];

    // *** DIAGRAM 10871 OF 15495 ***
    // Wavefunction(s) for diagram number 10871
    // (none)
    // Amplitude(s) for diagram number 10871
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[498], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];

    // *** DIAGRAM 10872 OF 15495 ***
    // Wavefunction(s) for diagram number 10872
    // (none)
    // Amplitude(s) for diagram number 10872
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[628], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];

    // *** DIAGRAM 10873 OF 15495 ***
    // Wavefunction(s) for diagram number 10873
    // (none)
    // Amplitude(s) for diagram number 10873
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[498], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10874 OF 15495 ***
    // Wavefunction(s) for diagram number 10874
    // (none)
    // Amplitude(s) for diagram number 10874
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[525], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];

    // *** DIAGRAM 10875 OF 15495 ***
    // Wavefunction(s) for diagram number 10875
    // (none)
    // Amplitude(s) for diagram number 10875
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[128], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 10876 OF 15495 ***
    // Wavefunction(s) for diagram number 10876
    // (none)
    // Amplitude(s) for diagram number 10876
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[697], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];

    // *** DIAGRAM 10877 OF 15495 ***
    // Wavefunction(s) for diagram number 10877
    // (none)
    // Amplitude(s) for diagram number 10877
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[274], w_fp[2], w_fp[130], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10878 OF 15495 ***
    // Wavefunction(s) for diagram number 10878
    // (none)
    // Amplitude(s) for diagram number 10878
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[525], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];

    // *** DIAGRAM 10879 OF 15495 ***
    // Wavefunction(s) for diagram number 10879
    // (none)
    // Amplitude(s) for diagram number 10879
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[128], w_fp[743], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];

    // *** DIAGRAM 10880 OF 15495 ***
    // Wavefunction(s) for diagram number 10880
    // (none)
    // Amplitude(s) for diagram number 10880
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[697], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[578] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[579] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[698] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[699] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10881 OF 15495 ***
    // Wavefunction(s) for diagram number 10881
    // (none)
    // Amplitude(s) for diagram number 10881
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[2], w_fp[743], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];

    // *** DIAGRAM 10882 OF 15495 ***
    // Wavefunction(s) for diagram number 10882
    // (none)
    // Amplitude(s) for diagram number 10882
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[628], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10883 OF 15495 ***
    // Wavefunction(s) for diagram number 10883
    // (none)
    // Amplitude(s) for diagram number 10883
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[128], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[578] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[128], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[578] += amp_sv[0];
    jamp_sv[579] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[587] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[597] -= amp_sv[0];
    jamp_sv[698] -= amp_sv[0];
    jamp_sv[699] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[707] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    jamp_sv[717] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[128], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[591] += amp_sv[0];
    jamp_sv[593] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 10884 OF 15495 ***
    // Wavefunction(s) for diagram number 10884
    // (none)
    // Amplitude(s) for diagram number 10884
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[2], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[2], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[30] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[188], w_fp[2], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];

    // *** DIAGRAM 10885 OF 15495 ***
    // Wavefunction(s) for diagram number 10885
    // (none)
    // Amplitude(s) for diagram number 10885
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[468], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[579], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[556], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[185] += amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 10886 OF 15495 ***
    // Wavefunction(s) for diagram number 10886
    // (none)
    // Amplitude(s) for diagram number 10886
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[468], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[579], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[556], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[46] -= amp_sv[0];
    jamp_sv[47] += amp_sv[0];
    jamp_sv[256] -= amp_sv[0];
    jamp_sv[257] += amp_sv[0];
    jamp_sv[262] += amp_sv[0];
    jamp_sv[263] -= amp_sv[0];
    jamp_sv[578] -= amp_sv[0];
    jamp_sv[579] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[698] += amp_sv[0];
    jamp_sv[699] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 10887 OF 15495 ***
    // Wavefunction(s) for diagram number 10887
    // (none)
    // Amplitude(s) for diagram number 10887
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[242], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[333] += amp_sv[0];
    jamp_sv[335] -= amp_sv[0];
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[713] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[47], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += amp_sv[0];
    jamp_sv[355] -= amp_sv[0];
    jamp_sv[357] -= amp_sv[0];
    jamp_sv[359] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[588] += amp_sv[0];
    jamp_sv[589] -= amp_sv[0];
    jamp_sv[591] -= amp_sv[0];
    jamp_sv[593] += amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[13], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[333] -= amp_sv[0];
    jamp_sv[335] += amp_sv[0];
    jamp_sv[546] -= amp_sv[0];
    jamp_sv[547] += amp_sv[0];
    jamp_sv[549] += amp_sv[0];
    jamp_sv[551] -= amp_sv[0];
    jamp_sv[666] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[671] -= amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    jamp_sv[713] += amp_sv[0];

    // *** DIAGRAM 10888 OF 15495 ***
    // Wavefunction(s) for diagram number 10888
    // (none)
    // Amplitude(s) for diagram number 10888
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[7] += amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[9] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 10889 OF 15495 ***
    // Wavefunction(s) for diagram number 10889
    // (none)
    // Amplitude(s) for diagram number 10889
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[593], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[581], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[185] -= amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[3], w_fp[494], w_fp[533], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[183] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 10890 OF 15495 ***
    // Wavefunction(s) for diagram number 10890
    // (none)
    // Amplitude(s) for diagram number 10890
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[593], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] -= amp_sv[0];
    jamp_sv[41] += amp_sv[0];
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[116] -= amp_sv[0];
    jamp_sv[117] += amp_sv[0];
    jamp_sv[330] += amp_sv[0];
    jamp_sv[331] -= amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[708] += amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[581], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[46] += amp_sv[0];
    jamp_sv[47] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[92] += amp_sv[0];
    jamp_sv[93] -= amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[354] -= amp_sv[0];
    jamp_sv[355] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[588] -= amp_sv[0];
    jamp_sv[589] += amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[2], w_fp[533], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[40] += amp_sv[0];
    jamp_sv[41] -= amp_sv[0];
    jamp_sv[82] -= amp_sv[0];
    jamp_sv[83] += amp_sv[0];
    jamp_sv[106] -= amp_sv[0];
    jamp_sv[107] += amp_sv[0];
    jamp_sv[116] += amp_sv[0];
    jamp_sv[117] -= amp_sv[0];
    jamp_sv[330] -= amp_sv[0];
    jamp_sv[331] += amp_sv[0];
    jamp_sv[546] += amp_sv[0];
    jamp_sv[547] -= amp_sv[0];
    jamp_sv[666] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[708] -= amp_sv[0];
    jamp_sv[709] += amp_sv[0];

    // *** DIAGRAM 10891 OF 15495 ***
    // Wavefunction(s) for diagram number 10891
    // (none)
    // Amplitude(s) for diagram number 10891
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[37], w_fp[494], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[40], w_fp[494], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[12], w_fp[494], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10892 OF 15495 ***
    // Wavefunction(s) for diagram number 10892
    // (none)
    // Amplitude(s) for diagram number 10892
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[242], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[47], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[354] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[355] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[588] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[589] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[302], w_fp[13], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[330] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[331] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[546] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[547] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[666] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[667] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[708] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10893 OF 15495 ***
    // Wavefunction(s) for diagram number 10893
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[0], w_fp[304], COUPs[0], 1.0, 0., 0., w_fp[662] );
    // Amplitude(s) for diagram number 10893
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[27], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[17] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[27], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[27], w_fp[5], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10894 OF 15495 ***
    // Wavefunction(s) for diagram number 10894
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[662], w_fp[5], COUPs[0], 1.0, 0., 0., w_fp[710] );
    // Amplitude(s) for diagram number 10894
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[7], w_fp[710], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[15] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[75] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[349] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[351] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[353] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10895 OF 15495 ***
    // Wavefunction(s) for diagram number 10895
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[662], w_fp[7], COUPs[0], 1.0, 0., 0., w_fp[711] );
    // Amplitude(s) for diagram number 10895
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[27], w_fp[5], w_fp[711], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[17] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[77] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[306] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[307] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[309] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[311] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[509] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10896 OF 15495 ***
    // Wavefunction(s) for diagram number 10896
    // (none)
    // Amplitude(s) for diagram number 10896
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[11], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[16] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[11], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[11], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10897 OF 15495 ***
    // Wavefunction(s) for diagram number 10897
    VVV1P0_1<W_ACCESS, CD_ACCESS>( w_fp[662], w_fp[4], COUPs[0], 1.0, 0., 0., w_fp[544] );
    // Amplitude(s) for diagram number 10897
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[11], w_fp[7], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[13] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[73] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[193] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[469] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[471] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[505] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10898 OF 15495 ***
    // Wavefunction(s) for diagram number 10898
    // (none)
    // Amplitude(s) for diagram number 10898
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[11], w_fp[4], w_fp[711], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[16] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[76] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[243] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[253] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[277] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[319] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[337] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[340] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[346] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[356] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[427] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[429] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[460] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[466] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[508] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10899 OF 15495 ***
    // Wavefunction(s) for diagram number 10899
    // (none)
    // Amplitude(s) for diagram number 10899
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[16], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[14] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] += cxtype( 0, 1 ) * amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[16], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[16], w_fp[4], w_fp[5], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[14] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[288] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[290] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[302] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[506] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[605] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[619] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[643] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[661] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[663] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[665] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10900 OF 15495 ***
    // Wavefunction(s) for diagram number 10900
    // (none)
    // Amplitude(s) for diagram number 10900
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[16], w_fp[5], w_fp[544], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[12] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[72] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[242] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[252] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[276] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[318] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[410] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[416] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[504] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[604] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[618] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[642] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[649] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[652] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[658] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[682] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[684] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[685] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[687] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[689] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[692] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[702] -= cxtype( 0, 1 ) * amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 544 );
    storeWf( wfs, w_cx, nevt, 662 );
    storeWf( wfs, w_cx, nevt, 710 );
    storeWf( wfs, w_cx, nevt, 711 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
