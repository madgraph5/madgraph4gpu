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
  diagramgroup65( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 10 );
    retrieveWf( wfs, w_cx, nevt, 12 );
    retrieveWf( wfs, w_cx, nevt, 17 );
    retrieveWf( wfs, w_cx, nevt, 18 );
    retrieveWf( wfs, w_cx, nevt, 19 );
    retrieveWf( wfs, w_cx, nevt, 20 );
    retrieveWf( wfs, w_cx, nevt, 21 );
    retrieveWf( wfs, w_cx, nevt, 22 );
    retrieveWf( wfs, w_cx, nevt, 26 );
    retrieveWf( wfs, w_cx, nevt, 37 );
    retrieveWf( wfs, w_cx, nevt, 40 );
    retrieveWf( wfs, w_cx, nevt, 45 );
    retrieveWf( wfs, w_cx, nevt, 54 );
    retrieveWf( wfs, w_cx, nevt, 61 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 71 );
    retrieveWf( wfs, w_cx, nevt, 80 );
    retrieveWf( wfs, w_cx, nevt, 83 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 104 );
    retrieveWf( wfs, w_cx, nevt, 105 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 132 );
    retrieveWf( wfs, w_cx, nevt, 159 );
    retrieveWf( wfs, w_cx, nevt, 160 );
    retrieveWf( wfs, w_cx, nevt, 167 );
    retrieveWf( wfs, w_cx, nevt, 169 );
    retrieveWf( wfs, w_cx, nevt, 173 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 177 );
    retrieveWf( wfs, w_cx, nevt, 178 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 187 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 189 );
    retrieveWf( wfs, w_cx, nevt, 190 );
    retrieveWf( wfs, w_cx, nevt, 191 );
    retrieveWf( wfs, w_cx, nevt, 193 );
    retrieveWf( wfs, w_cx, nevt, 194 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 201 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 203 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 205 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 207 );
    retrieveWf( wfs, w_cx, nevt, 208 );
    retrieveWf( wfs, w_cx, nevt, 209 );
    retrieveWf( wfs, w_cx, nevt, 210 );
    retrieveWf( wfs, w_cx, nevt, 212 );
    retrieveWf( wfs, w_cx, nevt, 213 );
    retrieveWf( wfs, w_cx, nevt, 216 );
    retrieveWf( wfs, w_cx, nevt, 238 );
    retrieveWf( wfs, w_cx, nevt, 248 );
    retrieveWf( wfs, w_cx, nevt, 252 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 264 );
    retrieveWf( wfs, w_cx, nevt, 291 );
    retrieveWf( wfs, w_cx, nevt, 305 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 363 );
    retrieveWf( wfs, w_cx, nevt, 367 );
    retrieveWf( wfs, w_cx, nevt, 368 );
    retrieveWf( wfs, w_cx, nevt, 369 );
    retrieveWf( wfs, w_cx, nevt, 371 );
    retrieveWf( wfs, w_cx, nevt, 372 );
    retrieveWf( wfs, w_cx, nevt, 373 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 476 );
    retrieveWf( wfs, w_cx, nevt, 477 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 486 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 529 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 540 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 579 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 593 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 604 );
    retrieveWf( wfs, w_cx, nevt, 606 );
    retrieveWf( wfs, w_cx, nevt, 607 );
    retrieveWf( wfs, w_cx, nevt, 609 );
    retrieveWf( wfs, w_cx, nevt, 610 );
    retrieveWf( wfs, w_cx, nevt, 611 );
    retrieveWf( wfs, w_cx, nevt, 613 );
    retrieveWf( wfs, w_cx, nevt, 614 );
    retrieveWf( wfs, w_cx, nevt, 616 );
    retrieveWf( wfs, w_cx, nevt, 617 );
    retrieveWf( wfs, w_cx, nevt, 628 );
    retrieveWf( wfs, w_cx, nevt, 639 );
    retrieveWf( wfs, w_cx, nevt, 640 );
    retrieveWf( wfs, w_cx, nevt, 641 );
    retrieveWf( wfs, w_cx, nevt, 642 );
    retrieveWf( wfs, w_cx, nevt, 643 );
    retrieveWf( wfs, w_cx, nevt, 644 );
    retrieveWf( wfs, w_cx, nevt, 676 );
    retrieveWf( wfs, w_cx, nevt, 677 );
    retrieveWf( wfs, w_cx, nevt, 678 );
    retrieveWf( wfs, w_cx, nevt, 679 );
    retrieveWf( wfs, w_cx, nevt, 681 );
    retrieveWf( wfs, w_cx, nevt, 682 );
    retrieveWf( wfs, w_cx, nevt, 683 );
    retrieveWf( wfs, w_cx, nevt, 684 );
    retrieveWf( wfs, w_cx, nevt, 685 );
    retrieveWf( wfs, w_cx, nevt, 700 );
    retrieveWf( wfs, w_cx, nevt, 742 );
    retrieveWf( wfs, w_cx, nevt, 744 );
#endif
#endif

    // *** DIAGRAM 12801 OF 15495 ***
    // Wavefunction(s) for diagram number 12801
    // (none)
    // Amplitude(s) for diagram number 12801
    FFV1_0( w_fp[196], w_fp[477], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];

    // *** DIAGRAM 12802 OF 15495 ***
    // Wavefunction(s) for diagram number 12802
    // (none)
    // Amplitude(s) for diagram number 12802
    FFV1_0( w_fp[355], w_fp[169], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 12803 OF 15495 ***
    // Wavefunction(s) for diagram number 12803
    // (none)
    // Amplitude(s) for diagram number 12803
    VVV1_0( w_fp[516], w_fp[1], w_fp[201], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12804 OF 15495 ***
    // Wavefunction(s) for diagram number 12804
    // (none)
    // Amplitude(s) for diagram number 12804
    FFV1_0( w_fp[129], w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12805 OF 15495 ***
    // Wavefunction(s) for diagram number 12805
    // (none)
    // Amplitude(s) for diagram number 12805
    FFV1_0( w_fp[355], w_fp[210], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[450] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12806 OF 15495 ***
    // Wavefunction(s) for diagram number 12806
    // (none)
    // Amplitude(s) for diagram number 12806
    VVV1_0( w_fp[0], w_fp[361], w_fp[201], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12807 OF 15495 ***
    // Wavefunction(s) for diagram number 12807
    // (none)
    // Amplitude(s) for diagram number 12807
    FFV1_0( w_fp[196], w_fp[169], w_fp[71], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[169], w_fp[132], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[196], w_fp[169], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[450] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[452] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[474] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[476] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12808 OF 15495 ***
    // Wavefunction(s) for diagram number 12808
    // (none)
    // Amplitude(s) for diagram number 12808
    FFV1_0( w_fp[255], w_fp[679], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] -= amp_sv[0];

    // *** DIAGRAM 12809 OF 15495 ***
    // Wavefunction(s) for diagram number 12809
    // (none)
    // Amplitude(s) for diagram number 12809
    FFV1_0( w_fp[255], w_fp[678], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[380] -= amp_sv[0];

    // *** DIAGRAM 12810 OF 15495 ***
    // Wavefunction(s) for diagram number 12810
    // (none)
    // Amplitude(s) for diagram number 12810
    FFV1_0( w_fp[202], w_fp[105], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] -= amp_sv[0];

    // *** DIAGRAM 12811 OF 15495 ***
    // Wavefunction(s) for diagram number 12811
    // (none)
    // Amplitude(s) for diagram number 12811
    FFV1_0( w_fp[202], w_fp[678], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= amp_sv[0];

    // *** DIAGRAM 12812 OF 15495 ***
    // Wavefunction(s) for diagram number 12812
    // (none)
    // Amplitude(s) for diagram number 12812
    FFV1_0( w_fp[176], w_fp[105], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];

    // *** DIAGRAM 12813 OF 15495 ***
    // Wavefunction(s) for diagram number 12813
    // (none)
    // Amplitude(s) for diagram number 12813
    FFV1_0( w_fp[176], w_fp[679], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[367] -= amp_sv[0];

    // *** DIAGRAM 12814 OF 15495 ***
    // Wavefunction(s) for diagram number 12814
    // (none)
    // Amplitude(s) for diagram number 12814
    FFV1_0( w_fp[595], w_fp[477], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[404] -= amp_sv[0];

    // *** DIAGRAM 12815 OF 15495 ***
    // Wavefunction(s) for diagram number 12815
    // (none)
    // Amplitude(s) for diagram number 12815
    FFV1_0( w_fp[578], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] -= amp_sv[0];

    // *** DIAGRAM 12816 OF 15495 ***
    // Wavefunction(s) for diagram number 12816
    // (none)
    // Amplitude(s) for diagram number 12816
    FFV1_0( w_fp[173], w_fp[190], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[428] -= amp_sv[0];

    // *** DIAGRAM 12817 OF 15495 ***
    // Wavefunction(s) for diagram number 12817
    // (none)
    // Amplitude(s) for diagram number 12817
    FFV1_0( w_fp[578], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[418] -= amp_sv[0];

    // *** DIAGRAM 12818 OF 15495 ***
    // Wavefunction(s) for diagram number 12818
    // (none)
    // Amplitude(s) for diagram number 12818
    FFV1_0( w_fp[173], w_fp[193], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[470] -= amp_sv[0];

    // *** DIAGRAM 12819 OF 15495 ***
    // Wavefunction(s) for diagram number 12819
    // (none)
    // Amplitude(s) for diagram number 12819
    FFV1_0( w_fp[595], w_fp[193], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[464] -= amp_sv[0];

    // *** DIAGRAM 12820 OF 15495 ***
    // Wavefunction(s) for diagram number 12820
    // (none)
    // Amplitude(s) for diagram number 12820
    FFV1_0( w_fp[202], w_fp[26], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] -= amp_sv[0];

    // *** DIAGRAM 12821 OF 15495 ***
    // Wavefunction(s) for diagram number 12821
    // (none)
    // Amplitude(s) for diagram number 12821
    FFV1_0( w_fp[545], w_fp[477], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[402] -= amp_sv[0];

    // *** DIAGRAM 12822 OF 15495 ***
    // Wavefunction(s) for diagram number 12822
    // (none)
    // Amplitude(s) for diagram number 12822
    FFV1_0( w_fp[176], w_fp[26], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] -= amp_sv[0];

    // *** DIAGRAM 12823 OF 15495 ***
    // Wavefunction(s) for diagram number 12823
    // (none)
    // Amplitude(s) for diagram number 12823
    FFV1_0( w_fp[555], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[391] -= amp_sv[0];

    // *** DIAGRAM 12824 OF 15495 ***
    // Wavefunction(s) for diagram number 12824
    // (none)
    // Amplitude(s) for diagram number 12824
    FFV1_0( w_fp[189], w_fp[190], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[426] -= amp_sv[0];

    // *** DIAGRAM 12825 OF 15495 ***
    // Wavefunction(s) for diagram number 12825
    // (none)
    // Amplitude(s) for diagram number 12825
    FFV1_0( w_fp[255], w_fp[683], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[412] -= amp_sv[0];

    // *** DIAGRAM 12826 OF 15495 ***
    // Wavefunction(s) for diagram number 12826
    // (none)
    // Amplitude(s) for diagram number 12826
    FFV1_0( w_fp[189], w_fp[193], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[468] -= amp_sv[0];

    // *** DIAGRAM 12827 OF 15495 ***
    // Wavefunction(s) for diagram number 12827
    // (none)
    // Amplitude(s) for diagram number 12827
    FFV1_0( w_fp[255], w_fp[682], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[458] -= amp_sv[0];

    // *** DIAGRAM 12828 OF 15495 ***
    // Wavefunction(s) for diagram number 12828
    // (none)
    // Amplitude(s) for diagram number 12828
    FFV1_0( w_fp[176], w_fp[683], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[409] -= amp_sv[0];

    // *** DIAGRAM 12829 OF 15495 ***
    // Wavefunction(s) for diagram number 12829
    // (none)
    // Amplitude(s) for diagram number 12829
    FFV1_0( w_fp[555], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[415] -= amp_sv[0];

    // *** DIAGRAM 12830 OF 15495 ***
    // Wavefunction(s) for diagram number 12830
    // (none)
    // Amplitude(s) for diagram number 12830
    FFV1_0( w_fp[202], w_fp[682], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= amp_sv[0];

    // *** DIAGRAM 12831 OF 15495 ***
    // Wavefunction(s) for diagram number 12831
    // (none)
    // Amplitude(s) for diagram number 12831
    FFV1_0( w_fp[545], w_fp[193], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] -= amp_sv[0];

    // *** DIAGRAM 12832 OF 15495 ***
    // Wavefunction(s) for diagram number 12832
    // (none)
    // Amplitude(s) for diagram number 12832
    FFV1_0( w_fp[255], w_fp[676], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12833 OF 15495 ***
    // Wavefunction(s) for diagram number 12833
    // (none)
    // Amplitude(s) for diagram number 12833
    FFV1_0( w_fp[174], w_fp[676], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];

    // *** DIAGRAM 12834 OF 15495 ***
    // Wavefunction(s) for diagram number 12834
    // (none)
    // Amplitude(s) for diagram number 12834
    FFV1_0( w_fp[83], w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12835 OF 15495 ***
    // Wavefunction(s) for diagram number 12835
    // (none)
    // Amplitude(s) for diagram number 12835
    FFV1_0( w_fp[594], w_fp[477], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12836 OF 15495 ***
    // Wavefunction(s) for diagram number 12836
    // (none)
    // Amplitude(s) for diagram number 12836
    FFV1_0( w_fp[594], w_fp[169], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 12837 OF 15495 ***
    // Wavefunction(s) for diagram number 12837
    // (none)
    // Amplitude(s) for diagram number 12837
    FFV1_0( w_fp[594], w_fp[209], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12838 OF 15495 ***
    // Wavefunction(s) for diagram number 12838
    // (none)
    // Amplitude(s) for diagram number 12838
    FFV1_0( w_fp[174], w_fp[477], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];

    // *** DIAGRAM 12839 OF 15495 ***
    // Wavefunction(s) for diagram number 12839
    // (none)
    // Amplitude(s) for diagram number 12839
    FFV1_0( w_fp[255], w_fp[169], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];

    // *** DIAGRAM 12840 OF 15495 ***
    // Wavefunction(s) for diagram number 12840
    // (none)
    // Amplitude(s) for diagram number 12840
    VVV1_0( w_fp[518], w_fp[1], w_fp[203], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12841 OF 15495 ***
    // Wavefunction(s) for diagram number 12841
    // (none)
    // Amplitude(s) for diagram number 12841
    FFV1_0( w_fp[83], w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12842 OF 15495 ***
    // Wavefunction(s) for diagram number 12842
    // (none)
    // Amplitude(s) for diagram number 12842
    FFV1_0( w_fp[255], w_fp[209], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[426] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12843 OF 15495 ***
    // Wavefunction(s) for diagram number 12843
    // (none)
    // Amplitude(s) for diagram number 12843
    VVV1_0( w_fp[0], w_fp[363], w_fp[203], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12844 OF 15495 ***
    // Wavefunction(s) for diagram number 12844
    // (none)
    // Amplitude(s) for diagram number 12844
    FFV1_0( w_fp[174], w_fp[169], w_fp[19], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[18], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[169], w_fp[17], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[426] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[428] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[468] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[470] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12845 OF 15495 ***
    // Wavefunction(s) for diagram number 12845
    // (none)
    // Amplitude(s) for diagram number 12845
    FFV1_0( w_fp[252], w_fp[679], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] -= amp_sv[0];

    // *** DIAGRAM 12846 OF 15495 ***
    // Wavefunction(s) for diagram number 12846
    // (none)
    // Amplitude(s) for diagram number 12846
    FFV1_0( w_fp[252], w_fp[677], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[374] -= amp_sv[0];

    // *** DIAGRAM 12847 OF 15495 ***
    // Wavefunction(s) for diagram number 12847
    // (none)
    // Amplitude(s) for diagram number 12847
    FFV1_0( w_fp[204], w_fp[105], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] -= amp_sv[0];

    // *** DIAGRAM 12848 OF 15495 ***
    // Wavefunction(s) for diagram number 12848
    // (none)
    // Amplitude(s) for diagram number 12848
    FFV1_0( w_fp[204], w_fp[677], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] -= amp_sv[0];

    // *** DIAGRAM 12849 OF 15495 ***
    // Wavefunction(s) for diagram number 12849
    // (none)
    // Amplitude(s) for diagram number 12849
    FFV1_0( w_fp[181], w_fp[105], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];

    // *** DIAGRAM 12850 OF 15495 ***
    // Wavefunction(s) for diagram number 12850
    // (none)
    // Amplitude(s) for diagram number 12850
    FFV1_0( w_fp[181], w_fp[679], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[366] -= amp_sv[0];

    // *** DIAGRAM 12851 OF 15495 ***
    // Wavefunction(s) for diagram number 12851
    // (none)
    // Amplitude(s) for diagram number 12851
    FFV1_0( w_fp[435], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[398] -= amp_sv[0];

    // *** DIAGRAM 12852 OF 15495 ***
    // Wavefunction(s) for diagram number 12852
    // (none)
    // Amplitude(s) for diagram number 12852
    FFV1_0( w_fp[495], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] -= amp_sv[0];

    // *** DIAGRAM 12853 OF 15495 ***
    // Wavefunction(s) for diagram number 12853
    // (none)
    // Amplitude(s) for diagram number 12853
    FFV1_0( w_fp[178], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[422] -= amp_sv[0];

    // *** DIAGRAM 12854 OF 15495 ***
    // Wavefunction(s) for diagram number 12854
    // (none)
    // Amplitude(s) for diagram number 12854
    FFV1_0( w_fp[495], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[416] -= amp_sv[0];

    // *** DIAGRAM 12855 OF 15495 ***
    // Wavefunction(s) for diagram number 12855
    // (none)
    // Amplitude(s) for diagram number 12855
    FFV1_0( w_fp[178], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[446] -= amp_sv[0];

    // *** DIAGRAM 12856 OF 15495 ***
    // Wavefunction(s) for diagram number 12856
    // (none)
    // Amplitude(s) for diagram number 12856
    FFV1_0( w_fp[435], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[440] -= amp_sv[0];

    // *** DIAGRAM 12857 OF 15495 ***
    // Wavefunction(s) for diagram number 12857
    // (none)
    // Amplitude(s) for diagram number 12857
    FFV1_0( w_fp[204], w_fp[26], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] -= amp_sv[0];

    // *** DIAGRAM 12858 OF 15495 ***
    // Wavefunction(s) for diagram number 12858
    // (none)
    // Amplitude(s) for diagram number 12858
    FFV1_0( w_fp[584], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[396] -= amp_sv[0];

    // *** DIAGRAM 12859 OF 15495 ***
    // Wavefunction(s) for diagram number 12859
    // (none)
    // Amplitude(s) for diagram number 12859
    FFV1_0( w_fp[181], w_fp[26], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= amp_sv[0];

    // *** DIAGRAM 12860 OF 15495 ***
    // Wavefunction(s) for diagram number 12860
    // (none)
    // Amplitude(s) for diagram number 12860
    FFV1_0( w_fp[583], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] -= amp_sv[0];

    // *** DIAGRAM 12861 OF 15495 ***
    // Wavefunction(s) for diagram number 12861
    // (none)
    // Amplitude(s) for diagram number 12861
    FFV1_0( w_fp[187], w_fp[190], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] -= amp_sv[0];

    // *** DIAGRAM 12862 OF 15495 ***
    // Wavefunction(s) for diagram number 12862
    // (none)
    // Amplitude(s) for diagram number 12862
    FFV1_0( w_fp[252], w_fp[683], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[410] -= amp_sv[0];

    // *** DIAGRAM 12863 OF 15495 ***
    // Wavefunction(s) for diagram number 12863
    // (none)
    // Amplitude(s) for diagram number 12863
    FFV1_0( w_fp[187], w_fp[191], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[444] -= amp_sv[0];

    // *** DIAGRAM 12864 OF 15495 ***
    // Wavefunction(s) for diagram number 12864
    // (none)
    // Amplitude(s) for diagram number 12864
    FFV1_0( w_fp[252], w_fp[681], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[434] -= amp_sv[0];

    // *** DIAGRAM 12865 OF 15495 ***
    // Wavefunction(s) for diagram number 12865
    // (none)
    // Amplitude(s) for diagram number 12865
    FFV1_0( w_fp[181], w_fp[683], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] -= amp_sv[0];

    // *** DIAGRAM 12866 OF 15495 ***
    // Wavefunction(s) for diagram number 12866
    // (none)
    // Amplitude(s) for diagram number 12866
    FFV1_0( w_fp[583], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[414] -= amp_sv[0];

    // *** DIAGRAM 12867 OF 15495 ***
    // Wavefunction(s) for diagram number 12867
    // (none)
    // Amplitude(s) for diagram number 12867
    FFV1_0( w_fp[204], w_fp[681], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= amp_sv[0];

    // *** DIAGRAM 12868 OF 15495 ***
    // Wavefunction(s) for diagram number 12868
    // (none)
    // Amplitude(s) for diagram number 12868
    FFV1_0( w_fp[584], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[438] -= amp_sv[0];

    // *** DIAGRAM 12869 OF 15495 ***
    // Wavefunction(s) for diagram number 12869
    // (none)
    // Amplitude(s) for diagram number 12869
    FFV1_0( w_fp[252], w_fp[676], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12870 OF 15495 ***
    // Wavefunction(s) for diagram number 12870
    // (none)
    // Amplitude(s) for diagram number 12870
    FFV1_0( w_fp[179], w_fp[676], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];

    // *** DIAGRAM 12871 OF 15495 ***
    // Wavefunction(s) for diagram number 12871
    // (none)
    // Amplitude(s) for diagram number 12871
    FFV1_0( w_fp[80], w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12872 OF 15495 ***
    // Wavefunction(s) for diagram number 12872
    // (none)
    // Amplitude(s) for diagram number 12872
    FFV1_0( w_fp[596], w_fp[477], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12873 OF 15495 ***
    // Wavefunction(s) for diagram number 12873
    // (none)
    // Amplitude(s) for diagram number 12873
    FFV1_0( w_fp[596], w_fp[169], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];

    // *** DIAGRAM 12874 OF 15495 ***
    // Wavefunction(s) for diagram number 12874
    // (none)
    // Amplitude(s) for diagram number 12874
    FFV1_0( w_fp[596], w_fp[207], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12875 OF 15495 ***
    // Wavefunction(s) for diagram number 12875
    // (none)
    // Amplitude(s) for diagram number 12875
    FFV1_0( w_fp[179], w_fp[477], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];

    // *** DIAGRAM 12876 OF 15495 ***
    // Wavefunction(s) for diagram number 12876
    // (none)
    // Amplitude(s) for diagram number 12876
    FFV1_0( w_fp[252], w_fp[169], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];

    // *** DIAGRAM 12877 OF 15495 ***
    // Wavefunction(s) for diagram number 12877
    // (none)
    // Amplitude(s) for diagram number 12877
    VVV1_0( w_fp[576], w_fp[1], w_fp[205], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12878 OF 15495 ***
    // Wavefunction(s) for diagram number 12878
    // (none)
    // Amplitude(s) for diagram number 12878
    FFV1_0( w_fp[80], w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12879 OF 15495 ***
    // Wavefunction(s) for diagram number 12879
    // (none)
    // Amplitude(s) for diagram number 12879
    FFV1_0( w_fp[252], w_fp[207], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[420] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12880 OF 15495 ***
    // Wavefunction(s) for diagram number 12880
    // (none)
    // Amplitude(s) for diagram number 12880
    VVV1_0( w_fp[0], w_fp[362], w_fp[205], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12881 OF 15495 ***
    // Wavefunction(s) for diagram number 12881
    // (none)
    // Amplitude(s) for diagram number 12881
    FFV1_0( w_fp[179], w_fp[169], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[169], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[169], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[420] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[422] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[444] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[446] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12882 OF 15495 ***
    // Wavefunction(s) for diagram number 12882
    // (none)
    // Amplitude(s) for diagram number 12882
    VVV1_0( w_fp[362], w_fp[7], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12883 OF 15495 ***
    // Wavefunction(s) for diagram number 12883
    // (none)
    // Amplitude(s) for diagram number 12883
    FFV1_0( w_fp[3], w_fp[678], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];

    // *** DIAGRAM 12884 OF 15495 ***
    // Wavefunction(s) for diagram number 12884
    // (none)
    // Amplitude(s) for diagram number 12884
    FFV1_0( w_fp[206], w_fp[105], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12885 OF 15495 ***
    // Wavefunction(s) for diagram number 12885
    // (none)
    // Amplitude(s) for diagram number 12885
    FFV1_0( w_fp[206], w_fp[678], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12886 OF 15495 ***
    // Wavefunction(s) for diagram number 12886
    // (none)
    // Amplitude(s) for diagram number 12886
    FFV1_0( w_fp[3], w_fp[105], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];

    // *** DIAGRAM 12887 OF 15495 ***
    // Wavefunction(s) for diagram number 12887
    // (none)
    // Amplitude(s) for diagram number 12887
    VVV1_0( w_fp[1], w_fp[88], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12888 OF 15495 ***
    // Wavefunction(s) for diagram number 12888
    // (none)
    // Amplitude(s) for diagram number 12888
    FFV1_0( w_fp[3], w_fp[676], w_fp[367], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[368], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[369], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[368] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[374] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12889 OF 15495 ***
    // Wavefunction(s) for diagram number 12889
    // (none)
    // Amplitude(s) for diagram number 12889
    FFV1_0( w_fp[607], w_fp[477], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];

    // *** DIAGRAM 12890 OF 15495 ***
    // Wavefunction(s) for diagram number 12890
    // (none)
    // Amplitude(s) for diagram number 12890
    FFV1_0( w_fp[3], w_fp[477], w_fp[606], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12891 OF 15495 ***
    // Wavefunction(s) for diagram number 12891
    // (none)
    // Amplitude(s) for diagram number 12891
    VVVV1_0( w_fp[576], w_fp[1], w_fp[194], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    VVVV3_0( w_fp[576], w_fp[1], w_fp[194], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    VVVV4_0( w_fp[576], w_fp[1], w_fp[194], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 12892 OF 15495 ***
    // Wavefunction(s) for diagram number 12892
    // (none)
    // Amplitude(s) for diagram number 12892
    VVV1_0( w_fp[194], w_fp[7], w_fp[604], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] += amp_sv[0];
    jamp_sv[374] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[420] -= amp_sv[0];
    jamp_sv[444] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 12893 OF 15495 ***
    // Wavefunction(s) for diagram number 12893
    // (none)
    // Amplitude(s) for diagram number 12893
    VVV1_0( w_fp[1], w_fp[194], w_fp[606], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[403] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[421] += amp_sv[0];
    jamp_sv[445] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 12894 OF 15495 ***
    // Wavefunction(s) for diagram number 12894
    // (none)
    // Amplitude(s) for diagram number 12894
    FFV1_0( w_fp[3], w_fp[193], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12895 OF 15495 ***
    // Wavefunction(s) for diagram number 12895
    // (none)
    // Amplitude(s) for diagram number 12895
    FFV1_0( w_fp[607], w_fp[193], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] += amp_sv[0];
    jamp_sv[463] -= amp_sv[0];
    jamp_sv[465] -= amp_sv[0];
    jamp_sv[467] += amp_sv[0];

    // *** DIAGRAM 12896 OF 15495 ***
    // Wavefunction(s) for diagram number 12896
    // (none)
    // Amplitude(s) for diagram number 12896
    FFV1_0( w_fp[206], w_fp[26], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12897 OF 15495 ***
    // Wavefunction(s) for diagram number 12897
    // (none)
    // Amplitude(s) for diagram number 12897
    FFV1_0( w_fp[610], w_fp[477], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[402] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12898 OF 15495 ***
    // Wavefunction(s) for diagram number 12898
    // (none)
    // Amplitude(s) for diagram number 12898
    FFV1_0( w_fp[3], w_fp[26], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 12899 OF 15495 ***
    // Wavefunction(s) for diagram number 12899
    // (none)
    // Amplitude(s) for diagram number 12899
    FFV1_0( w_fp[3], w_fp[477], w_fp[609], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12900 OF 15495 ***
    // Wavefunction(s) for diagram number 12900
    // (none)
    // Amplitude(s) for diagram number 12900
    VVVV1_0( w_fp[0], w_fp[362], w_fp[194], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[362], w_fp[194], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[362], w_fp[194], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];

    // *** DIAGRAM 12901 OF 15495 ***
    // Wavefunction(s) for diagram number 12901
    // (none)
    // Amplitude(s) for diagram number 12901
    VVV1_0( w_fp[194], w_fp[7], w_fp[64], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12902 OF 15495 ***
    // Wavefunction(s) for diagram number 12902
    // (none)
    // Amplitude(s) for diagram number 12902
    VVV1_0( w_fp[362], w_fp[7], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12903 OF 15495 ***
    // Wavefunction(s) for diagram number 12903
    // (none)
    // Amplitude(s) for diagram number 12903
    FFV1_0( w_fp[3], w_fp[193], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12904 OF 15495 ***
    // Wavefunction(s) for diagram number 12904
    // (none)
    // Amplitude(s) for diagram number 12904
    FFV1_0( w_fp[3], w_fp[682], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];

    // *** DIAGRAM 12905 OF 15495 ***
    // Wavefunction(s) for diagram number 12905
    // (none)
    // Amplitude(s) for diagram number 12905
    VVVV1_0( w_fp[0], w_fp[1], w_fp[194], w_fp[88], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[194], w_fp[88], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[194], w_fp[88], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 12906 OF 15495 ***
    // Wavefunction(s) for diagram number 12906
    // (none)
    // Amplitude(s) for diagram number 12906
    VVV1_0( w_fp[1], w_fp[88], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12907 OF 15495 ***
    // Wavefunction(s) for diagram number 12907
    // (none)
    // Amplitude(s) for diagram number 12907
    VVV1_0( w_fp[1], w_fp[194], w_fp[609], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 12908 OF 15495 ***
    // Wavefunction(s) for diagram number 12908
    // (none)
    // Amplitude(s) for diagram number 12908
    FFV1_0( w_fp[206], w_fp[682], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12909 OF 15495 ***
    // Wavefunction(s) for diagram number 12909
    // (none)
    // Amplitude(s) for diagram number 12909
    FFV1_0( w_fp[610], w_fp[193], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[462] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12910 OF 15495 ***
    // Wavefunction(s) for diagram number 12910
    // (none)
    // Amplitude(s) for diagram number 12910
    VVV1_0( w_fp[22], w_fp[194], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[422] += amp_sv[0];
    jamp_sv[446] -= amp_sv[0];
    jamp_sv[456] -= amp_sv[0];
    jamp_sv[457] += amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVV1_0( w_fp[21], w_fp[194], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[459] += amp_sv[0];
    jamp_sv[461] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    VVV1_0( w_fp[20], w_fp[194], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[420] += amp_sv[0];
    jamp_sv[422] -= amp_sv[0];
    jamp_sv[444] -= amp_sv[0];
    jamp_sv[446] += amp_sv[0];
    jamp_sv[456] += amp_sv[0];
    jamp_sv[457] -= amp_sv[0];
    jamp_sv[462] -= amp_sv[0];
    jamp_sv[463] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 12911 OF 15495 ***
    // Wavefunction(s) for diagram number 12911
    // (none)
    // Amplitude(s) for diagram number 12911
    FFV1_0( w_fp[3], w_fp[193], w_fp[22], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[193], w_fp[21], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[459] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[461] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[465] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[467] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[193], w_fp[20], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[456] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[457] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[462] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[463] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[472] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[473] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[478] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[479] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12912 OF 15495 ***
    // Wavefunction(s) for diagram number 12912
    // (none)
    // Amplitude(s) for diagram number 12912
    FFV1_0( w_fp[3], w_fp[477], w_fp[550], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[403] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12913 OF 15495 ***
    // Wavefunction(s) for diagram number 12913
    // (none)
    // Amplitude(s) for diagram number 12913
    VVV1_0( w_fp[550], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    VVV1_0( w_fp[480], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    VVV1_0( w_fp[558], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[403] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[421] -= amp_sv[0];
    jamp_sv[445] += amp_sv[0];
    jamp_sv[459] -= amp_sv[0];
    jamp_sv[461] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 12914 OF 15495 ***
    // Wavefunction(s) for diagram number 12914
    // (none)
    // Amplitude(s) for diagram number 12914
    VVV1_0( w_fp[0], w_fp[367], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[368], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[369], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[368] -= amp_sv[0];
    jamp_sv[374] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[423] += amp_sv[0];
    jamp_sv[447] -= amp_sv[0];
    jamp_sv[465] += amp_sv[0];
    jamp_sv[467] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12915 OF 15495 ***
    // Wavefunction(s) for diagram number 12915
    // (none)
    // Amplitude(s) for diagram number 12915
    VVV1_0( w_fp[363], w_fp[6], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12916 OF 15495 ***
    // Wavefunction(s) for diagram number 12916
    // (none)
    // Amplitude(s) for diagram number 12916
    FFV1_0( w_fp[3], w_fp[677], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];

    // *** DIAGRAM 12917 OF 15495 ***
    // Wavefunction(s) for diagram number 12917
    // (none)
    // Amplitude(s) for diagram number 12917
    FFV1_0( w_fp[208], w_fp[105], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12918 OF 15495 ***
    // Wavefunction(s) for diagram number 12918
    // (none)
    // Amplitude(s) for diagram number 12918
    FFV1_0( w_fp[208], w_fp[677], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12919 OF 15495 ***
    // Wavefunction(s) for diagram number 12919
    // (none)
    // Amplitude(s) for diagram number 12919
    FFV1_0( w_fp[3], w_fp[105], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];

    // *** DIAGRAM 12920 OF 15495 ***
    // Wavefunction(s) for diagram number 12920
    // (none)
    // Amplitude(s) for diagram number 12920
    VVV1_0( w_fp[1], w_fp[104], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12921 OF 15495 ***
    // Wavefunction(s) for diagram number 12921
    // (none)
    // Amplitude(s) for diagram number 12921
    FFV1_0( w_fp[3], w_fp[676], w_fp[371], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[372], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[373], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[372] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12922 OF 15495 ***
    // Wavefunction(s) for diagram number 12922
    // (none)
    // Amplitude(s) for diagram number 12922
    FFV1_0( w_fp[614], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];

    // *** DIAGRAM 12923 OF 15495 ***
    // Wavefunction(s) for diagram number 12923
    // (none)
    // Amplitude(s) for diagram number 12923
    FFV1_0( w_fp[3], w_fp[477], w_fp[613], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12924 OF 15495 ***
    // Wavefunction(s) for diagram number 12924
    // (none)
    // Amplitude(s) for diagram number 12924
    VVVV1_0( w_fp[518], w_fp[1], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    VVVV3_0( w_fp[518], w_fp[1], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[468] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];
    VVVV4_0( w_fp[518], w_fp[1], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];

    // *** DIAGRAM 12925 OF 15495 ***
    // Wavefunction(s) for diagram number 12925
    // (none)
    // Amplitude(s) for diagram number 12925
    VVV1_0( w_fp[194], w_fp[6], w_fp[611], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[426] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[468] += amp_sv[0];

    // *** DIAGRAM 12926 OF 15495 ***
    // Wavefunction(s) for diagram number 12926
    // (none)
    // Amplitude(s) for diagram number 12926
    VVV1_0( w_fp[1], w_fp[194], w_fp[613], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[397] += amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[427] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[469] -= amp_sv[0];

    // *** DIAGRAM 12927 OF 15495 ***
    // Wavefunction(s) for diagram number 12927
    // (none)
    // Amplitude(s) for diagram number 12927
    FFV1_0( w_fp[3], w_fp[191], w_fp[611], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12928 OF 15495 ***
    // Wavefunction(s) for diagram number 12928
    // (none)
    // Amplitude(s) for diagram number 12928
    FFV1_0( w_fp[614], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[438] += amp_sv[0];
    jamp_sv[439] -= amp_sv[0];
    jamp_sv[441] -= amp_sv[0];
    jamp_sv[443] += amp_sv[0];

    // *** DIAGRAM 12929 OF 15495 ***
    // Wavefunction(s) for diagram number 12929
    // (none)
    // Amplitude(s) for diagram number 12929
    FFV1_0( w_fp[208], w_fp[26], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12930 OF 15495 ***
    // Wavefunction(s) for diagram number 12930
    // (none)
    // Amplitude(s) for diagram number 12930
    FFV1_0( w_fp[617], w_fp[477], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[396] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12931 OF 15495 ***
    // Wavefunction(s) for diagram number 12931
    // (none)
    // Amplitude(s) for diagram number 12931
    FFV1_0( w_fp[3], w_fp[26], w_fp[104], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];

    // *** DIAGRAM 12932 OF 15495 ***
    // Wavefunction(s) for diagram number 12932
    // (none)
    // Amplitude(s) for diagram number 12932
    FFV1_0( w_fp[3], w_fp[477], w_fp[616], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12933 OF 15495 ***
    // Wavefunction(s) for diagram number 12933
    // (none)
    // Amplitude(s) for diagram number 12933
    VVVV1_0( w_fp[0], w_fp[363], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[363], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[363], w_fp[194], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[470] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];

    // *** DIAGRAM 12934 OF 15495 ***
    // Wavefunction(s) for diagram number 12934
    // (none)
    // Amplitude(s) for diagram number 12934
    VVV1_0( w_fp[194], w_fp[6], w_fp[61], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];

    // *** DIAGRAM 12935 OF 15495 ***
    // Wavefunction(s) for diagram number 12935
    // (none)
    // Amplitude(s) for diagram number 12935
    VVV1_0( w_fp[363], w_fp[6], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];

    // *** DIAGRAM 12936 OF 15495 ***
    // Wavefunction(s) for diagram number 12936
    // (none)
    // Amplitude(s) for diagram number 12936
    FFV1_0( w_fp[3], w_fp[191], w_fp[61], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12937 OF 15495 ***
    // Wavefunction(s) for diagram number 12937
    // (none)
    // Amplitude(s) for diagram number 12937
    FFV1_0( w_fp[3], w_fp[681], w_fp[363], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];

    // *** DIAGRAM 12938 OF 15495 ***
    // Wavefunction(s) for diagram number 12938
    // (none)
    // Amplitude(s) for diagram number 12938
    VVVV1_0( w_fp[0], w_fp[1], w_fp[194], w_fp[104], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[194], w_fp[104], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[194], w_fp[104], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];

    // *** DIAGRAM 12939 OF 15495 ***
    // Wavefunction(s) for diagram number 12939
    // (none)
    // Amplitude(s) for diagram number 12939
    VVV1_0( w_fp[1], w_fp[104], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];

    // *** DIAGRAM 12940 OF 15495 ***
    // Wavefunction(s) for diagram number 12940
    // (none)
    // Amplitude(s) for diagram number 12940
    VVV1_0( w_fp[1], w_fp[194], w_fp[616], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];

    // *** DIAGRAM 12941 OF 15495 ***
    // Wavefunction(s) for diagram number 12941
    // (none)
    // Amplitude(s) for diagram number 12941
    FFV1_0( w_fp[208], w_fp[681], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12942 OF 15495 ***
    // Wavefunction(s) for diagram number 12942
    // (none)
    // Amplitude(s) for diagram number 12942
    FFV1_0( w_fp[617], w_fp[191], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[438] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12943 OF 15495 ***
    // Wavefunction(s) for diagram number 12943
    // (none)
    // Amplitude(s) for diagram number 12943
    VVV1_0( w_fp[19], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[428] += amp_sv[0];
    jamp_sv[432] -= amp_sv[0];
    jamp_sv[433] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[470] -= amp_sv[0];
    VVV1_0( w_fp[18], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[435] += amp_sv[0];
    jamp_sv[437] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    VVV1_0( w_fp[17], w_fp[194], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[426] += amp_sv[0];
    jamp_sv[428] -= amp_sv[0];
    jamp_sv[432] += amp_sv[0];
    jamp_sv[433] -= amp_sv[0];
    jamp_sv[438] -= amp_sv[0];
    jamp_sv[439] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[468] -= amp_sv[0];
    jamp_sv[470] += amp_sv[0];

    // *** DIAGRAM 12944 OF 15495 ***
    // Wavefunction(s) for diagram number 12944
    // (none)
    // Amplitude(s) for diagram number 12944
    FFV1_0( w_fp[3], w_fp[191], w_fp[19], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[191], w_fp[18], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[435] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[437] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[441] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[443] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[191], w_fp[17], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[432] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[433] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[438] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[439] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[448] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[449] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[454] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[455] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12945 OF 15495 ***
    // Wavefunction(s) for diagram number 12945
    // (none)
    // Amplitude(s) for diagram number 12945
    FFV1_0( w_fp[3], w_fp[477], w_fp[486], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[547], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[529], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[397] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12946 OF 15495 ***
    // Wavefunction(s) for diagram number 12946
    // (none)
    // Amplitude(s) for diagram number 12946
    VVV1_0( w_fp[486], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    VVV1_0( w_fp[547], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[469] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    VVV1_0( w_fp[529], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[397] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[427] -= amp_sv[0];
    jamp_sv[435] -= amp_sv[0];
    jamp_sv[437] += amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[469] += amp_sv[0];

    // *** DIAGRAM 12947 OF 15495 ***
    // Wavefunction(s) for diagram number 12947
    // (none)
    // Amplitude(s) for diagram number 12947
    VVV1_0( w_fp[0], w_fp[371], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[372], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[471] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[373], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[370] -= amp_sv[0];
    jamp_sv[372] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[429] += amp_sv[0];
    jamp_sv[441] += amp_sv[0];
    jamp_sv[443] -= amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[471] -= amp_sv[0];

    // *** DIAGRAM 12948 OF 15495 ***
    // Wavefunction(s) for diagram number 12948
    // (none)
    // Amplitude(s) for diagram number 12948
    VVV1_0( w_fp[361], w_fp[4], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12949 OF 15495 ***
    // Wavefunction(s) for diagram number 12949
    // (none)
    // Amplitude(s) for diagram number 12949
    FFV1_0( w_fp[3], w_fp[679], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[366] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];

    // *** DIAGRAM 12950 OF 15495 ***
    // Wavefunction(s) for diagram number 12950
    // (none)
    // Amplitude(s) for diagram number 12950
    FFV1_0( w_fp[188], w_fp[105], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12951 OF 15495 ***
    // Wavefunction(s) for diagram number 12951
    // (none)
    // Amplitude(s) for diagram number 12951
    FFV1_0( w_fp[188], w_fp[679], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[366] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12952 OF 15495 ***
    // Wavefunction(s) for diagram number 12952
    // (none)
    // Amplitude(s) for diagram number 12952
    FFV1_0( w_fp[3], w_fp[105], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];

    // *** DIAGRAM 12953 OF 15495 ***
    // Wavefunction(s) for diagram number 12953
    // (none)
    // Amplitude(s) for diagram number 12953
    VVV1_0( w_fp[1], w_fp[130], w_fp[742], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12954 OF 15495 ***
    // Wavefunction(s) for diagram number 12954
    // (none)
    // Amplitude(s) for diagram number 12954
    FFV1_0( w_fp[3], w_fp[676], w_fp[291], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[248], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[264], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[366] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12955 OF 15495 ***
    // Wavefunction(s) for diagram number 12955
    // (none)
    // Amplitude(s) for diagram number 12955
    FFV1_0( w_fp[442], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];

    // *** DIAGRAM 12956 OF 15495 ***
    // Wavefunction(s) for diagram number 12956
    // (none)
    // Amplitude(s) for diagram number 12956
    FFV1_0( w_fp[3], w_fp[477], w_fp[546], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12957 OF 15495 ***
    // Wavefunction(s) for diagram number 12957
    // (none)
    // Amplitude(s) for diagram number 12957
    VVVV1_0( w_fp[516], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    VVVV3_0( w_fp[516], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[474] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];
    VVVV4_0( w_fp[516], w_fp[1], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 12958 OF 15495 ***
    // Wavefunction(s) for diagram number 12958
    // (none)
    // Amplitude(s) for diagram number 12958
    VVV1_0( w_fp[194], w_fp[4], w_fp[700], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[450] -= amp_sv[0];
    jamp_sv[474] += amp_sv[0];

    // *** DIAGRAM 12959 OF 15495 ***
    // Wavefunction(s) for diagram number 12959
    // (none)
    // Amplitude(s) for diagram number 12959
    VVV1_0( w_fp[1], w_fp[194], w_fp[546], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[451] += amp_sv[0];
    jamp_sv[475] -= amp_sv[0];

    // *** DIAGRAM 12960 OF 15495 ***
    // Wavefunction(s) for diagram number 12960
    // (none)
    // Amplitude(s) for diagram number 12960
    FFV1_0( w_fp[3], w_fp[190], w_fp[700], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12961 OF 15495 ***
    // Wavefunction(s) for diagram number 12961
    // (none)
    // Amplitude(s) for diagram number 12961
    FFV1_0( w_fp[442], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[414] += amp_sv[0];
    jamp_sv[415] -= amp_sv[0];
    jamp_sv[417] -= amp_sv[0];
    jamp_sv[419] += amp_sv[0];

    // *** DIAGRAM 12962 OF 15495 ***
    // Wavefunction(s) for diagram number 12962
    // (none)
    // Amplitude(s) for diagram number 12962
    FFV1_0( w_fp[188], w_fp[26], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12963 OF 15495 ***
    // Wavefunction(s) for diagram number 12963
    // (none)
    // Amplitude(s) for diagram number 12963
    FFV1_0( w_fp[628], w_fp[477], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[390] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12964 OF 15495 ***
    // Wavefunction(s) for diagram number 12964
    // (none)
    // Amplitude(s) for diagram number 12964
    FFV1_0( w_fp[3], w_fp[26], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 12965 OF 15495 ***
    // Wavefunction(s) for diagram number 12965
    // (none)
    // Amplitude(s) for diagram number 12965
    FFV1_0( w_fp[3], w_fp[477], w_fp[525], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12966 OF 15495 ***
    // Wavefunction(s) for diagram number 12966
    // (none)
    // Amplitude(s) for diagram number 12966
    VVVV1_0( w_fp[0], w_fp[361], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[361], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[361], w_fp[194], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[476] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];

    // *** DIAGRAM 12967 OF 15495 ***
    // Wavefunction(s) for diagram number 12967
    // (none)
    // Amplitude(s) for diagram number 12967
    VVV1_0( w_fp[194], w_fp[4], w_fp[238], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];

    // *** DIAGRAM 12968 OF 15495 ***
    // Wavefunction(s) for diagram number 12968
    // (none)
    // Amplitude(s) for diagram number 12968
    VVV1_0( w_fp[361], w_fp[4], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];

    // *** DIAGRAM 12969 OF 15495 ***
    // Wavefunction(s) for diagram number 12969
    // (none)
    // Amplitude(s) for diagram number 12969
    FFV1_0( w_fp[3], w_fp[190], w_fp[238], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12970 OF 15495 ***
    // Wavefunction(s) for diagram number 12970
    // (none)
    // Amplitude(s) for diagram number 12970
    FFV1_0( w_fp[3], w_fp[683], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];

    // *** DIAGRAM 12971 OF 15495 ***
    // Wavefunction(s) for diagram number 12971
    // (none)
    // Amplitude(s) for diagram number 12971
    VVVV1_0( w_fp[0], w_fp[1], w_fp[194], w_fp[130], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[1], w_fp[194], w_fp[130], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[1], w_fp[194], w_fp[130], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 12972 OF 15495 ***
    // Wavefunction(s) for diagram number 12972
    // (none)
    // Amplitude(s) for diagram number 12972
    VVV1_0( w_fp[1], w_fp[130], w_fp[744], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12973 OF 15495 ***
    // Wavefunction(s) for diagram number 12973
    // (none)
    // Amplitude(s) for diagram number 12973
    VVV1_0( w_fp[1], w_fp[194], w_fp[525], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];

    // *** DIAGRAM 12974 OF 15495 ***
    // Wavefunction(s) for diagram number 12974
    // (none)
    // Amplitude(s) for diagram number 12974
    FFV1_0( w_fp[188], w_fp[683], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12975 OF 15495 ***
    // Wavefunction(s) for diagram number 12975
    // (none)
    // Amplitude(s) for diagram number 12975
    FFV1_0( w_fp[628], w_fp[190], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[414] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12976 OF 15495 ***
    // Wavefunction(s) for diagram number 12976
    // (none)
    // Amplitude(s) for diagram number 12976
    VVV1_0( w_fp[71], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[408] -= amp_sv[0];
    jamp_sv[409] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[452] += amp_sv[0];
    jamp_sv[476] -= amp_sv[0];
    VVV1_0( w_fp[132], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[411] += amp_sv[0];
    jamp_sv[413] -= amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    VVV1_0( w_fp[65], w_fp[194], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[408] += amp_sv[0];
    jamp_sv[409] -= amp_sv[0];
    jamp_sv[414] -= amp_sv[0];
    jamp_sv[415] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[450] += amp_sv[0];
    jamp_sv[452] -= amp_sv[0];
    jamp_sv[474] -= amp_sv[0];
    jamp_sv[476] += amp_sv[0];

    // *** DIAGRAM 12977 OF 15495 ***
    // Wavefunction(s) for diagram number 12977
    // (none)
    // Amplitude(s) for diagram number 12977
    FFV1_0( w_fp[3], w_fp[190], w_fp[71], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[190], w_fp[132], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[411] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[413] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[417] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[419] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[190], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[408] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[409] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[414] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[415] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[424] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[425] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[430] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[431] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12978 OF 15495 ***
    // Wavefunction(s) for diagram number 12978
    // (none)
    // Amplitude(s) for diagram number 12978
    FFV1_0( w_fp[3], w_fp[477], w_fp[468], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12979 OF 15495 ***
    // Wavefunction(s) for diagram number 12979
    // (none)
    // Amplitude(s) for diagram number 12979
    VVV1_0( w_fp[468], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    VVV1_0( w_fp[579], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    VVV1_0( w_fp[556], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[411] -= amp_sv[0];
    jamp_sv[413] += amp_sv[0];
    jamp_sv[451] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[475] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 12980 OF 15495 ***
    // Wavefunction(s) for diagram number 12980
    // (none)
    // Amplitude(s) for diagram number 12980
    VVV1_0( w_fp[0], w_fp[291], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[248], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[264], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[366] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[417] += amp_sv[0];
    jamp_sv[419] -= amp_sv[0];
    jamp_sv[453] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[477] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12981 OF 15495 ***
    // Wavefunction(s) for diagram number 12981
    // (none)
    // Amplitude(s) for diagram number 12981
    FFV1_0( w_fp[3], w_fp[676], w_fp[177], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[305], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[371] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[377] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[676], w_fp[10], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[362] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[369] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[375] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[381] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[383] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12982 OF 15495 ***
    // Wavefunction(s) for diagram number 12982
    // (none)
    // Amplitude(s) for diagram number 12982
    FFV1_0( w_fp[37], w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    FFV1_0( w_fp[40], w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    FFV1_0( w_fp[12], w_fp[676], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];

    // *** DIAGRAM 12983 OF 15495 ***
    // Wavefunction(s) for diagram number 12983
    // (none)
    // Amplitude(s) for diagram number 12983
    FFV1_0( w_fp[3], w_fp[477], w_fp[593], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[395] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[401] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[477], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[393] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[399] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[405] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[407] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 12984 OF 15495 ***
    // Wavefunction(s) for diagram number 12984
    // (none)
    // Amplitude(s) for diagram number 12984
    VVV1_0( w_fp[593], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[371] += amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    VVV1_0( w_fp[581], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[377] += amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[395] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[401] += amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[430] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[454] -= amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    VVV1_0( w_fp[533], w_fp[1], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[375] -= amp_sv[0];
    jamp_sv[381] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[399] -= amp_sv[0];
    jamp_sv[405] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[448] += amp_sv[0];
    jamp_sv[472] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];

    // *** DIAGRAM 12985 OF 15495 ***
    // Wavefunction(s) for diagram number 12985
    // (none)
    // Amplitude(s) for diagram number 12985
    FFV1_0( w_fp[37], w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    FFV1_0( w_fp[40], w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    FFV1_0( w_fp[12], w_fp[477], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];

    // *** DIAGRAM 12986 OF 15495 ***
    // Wavefunction(s) for diagram number 12986
    // (none)
    // Amplitude(s) for diagram number 12986
    VVV1_0( w_fp[0], w_fp[177], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[305], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    VVV1_0( w_fp[0], w_fp[10], w_fp[194], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];

    // *** DIAGRAM 12987 OF 15495 ***
    // Wavefunction(s) for diagram number 12987
    // (none)
    // Amplitude(s) for diagram number 12987
    FFV1_0( w_fp[3], w_fp[169], w_fp[167], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[54], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] += amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[383] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[393] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[407] += amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[159], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[424] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[478] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[644], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[431] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[455] += amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[643], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[371] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[377] -= amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[395] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[401] -= amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[642], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[430] += amp_sv[0];
    jamp_sv[431] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[454] += amp_sv[0];
    jamp_sv[455] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[641], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] += amp_sv[0];
    jamp_sv[362] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[425] += amp_sv[0];
    jamp_sv[449] -= amp_sv[0];
    jamp_sv[473] -= amp_sv[0];
    jamp_sv[479] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[640], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[369] -= amp_sv[0];
    jamp_sv[375] += amp_sv[0];
    jamp_sv[381] += amp_sv[0];
    jamp_sv[383] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[393] -= amp_sv[0];
    jamp_sv[399] += amp_sv[0];
    jamp_sv[405] += amp_sv[0];
    jamp_sv[407] -= amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    FFV1_0( w_fp[3], w_fp[169], w_fp[639], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[360] -= amp_sv[0];
    jamp_sv[362] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[424] += amp_sv[0];
    jamp_sv[425] -= amp_sv[0];
    jamp_sv[448] -= amp_sv[0];
    jamp_sv[449] += amp_sv[0];
    jamp_sv[472] -= amp_sv[0];
    jamp_sv[473] += amp_sv[0];
    jamp_sv[478] += amp_sv[0];
    jamp_sv[479] -= amp_sv[0];

    // *** DIAGRAM 12988 OF 15495 ***
    // Wavefunction(s) for diagram number 12988
    // (none)
    // Amplitude(s) for diagram number 12988
    FFV1_0( w_fp[355], w_fp[684], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[496] -= amp_sv[0];

    // *** DIAGRAM 12989 OF 15495 ***
    // Wavefunction(s) for diagram number 12989
    // (none)
    // Amplitude(s) for diagram number 12989
    FFV1_0( w_fp[355], w_fp[685], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[502] -= amp_sv[0];

    // *** DIAGRAM 12990 OF 15495 ***
    // Wavefunction(s) for diagram number 12990
    FFV1_1( w_fp[540], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[12] );
    // Amplitude(s) for diagram number 12990
    FFV1_0( w_fp[216], w_fp[12], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[485] -= amp_sv[0];

    // *** DIAGRAM 12991 OF 15495 ***
    // Wavefunction(s) for diagram number 12991
    // (none)
    // Amplitude(s) for diagram number 12991
    FFV1_0( w_fp[216], w_fp[685], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[499] -= amp_sv[0];

    // *** DIAGRAM 12992 OF 15495 ***
    // Wavefunction(s) for diagram number 12992
    // (none)
    // Amplitude(s) for diagram number 12992
    FFV1_0( w_fp[199], w_fp[12], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[483] -= amp_sv[0];

    // *** DIAGRAM 12993 OF 15495 ***
    // Wavefunction(s) for diagram number 12993
    // (none)
    // Amplitude(s) for diagram number 12993
    FFV1_0( w_fp[199], w_fp[684], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[493] -= amp_sv[0];

    // *** DIAGRAM 12994 OF 15495 ***
    // Wavefunction(s) for diagram number 12994
    // (none)
    // Amplitude(s) for diagram number 12994
    FFV1_0( w_fp[586], w_fp[476], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[526] -= amp_sv[0];

    // *** DIAGRAM 12995 OF 15495 ***
    // Wavefunction(s) for diagram number 12995
    // (none)
    // Amplitude(s) for diagram number 12995
    FFV1_0( w_fp[45], w_fp[476], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[520] -= amp_sv[0];

    // *** DIAGRAM 12996 OF 15495 ***
    // Wavefunction(s) for diagram number 12996
    // (none)
    // Amplitude(s) for diagram number 12996
    FFV1_0( w_fp[160], w_fp[212], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[572] -= amp_sv[0];

    // *** DIAGRAM 12997 OF 15495 ***
    // Wavefunction(s) for diagram number 12997
    // (none)
    // Amplitude(s) for diagram number 12997
    FFV1_0( w_fp[45], w_fp[212], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[562] -= amp_sv[0];

    // *** DIAGRAM 12998 OF 15495 ***
    // Wavefunction(s) for diagram number 12998
    // (none)
    // Amplitude(s) for diagram number 12998
    FFV1_0( w_fp[160], w_fp[213], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[596] -= amp_sv[0];

    // *** DIAGRAM 12999 OF 15495 ***
    // Wavefunction(s) for diagram number 12999
    // (none)
    // Amplitude(s) for diagram number 12999
    FFV1_0( w_fp[586], w_fp[213], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[586] -= amp_sv[0];

    // *** DIAGRAM 13000 OF 15495 ***
    // Wavefunction(s) for diagram number 13000
    FFV1_1( w_fp[476], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[40] );
    // Amplitude(s) for diagram number 13000
    FFV1_0( w_fp[216], w_fp[40], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[509] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 12 );
    storeWf( wfs, w_cx, nevt, 40 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
