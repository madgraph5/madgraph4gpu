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
  diagramgroup60( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 198 );
    retrieveWf( wfs, w_cx, nevt, 199 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 215 );
    retrieveWf( wfs, w_cx, nevt, 225 );
    retrieveWf( wfs, w_cx, nevt, 228 );
    retrieveWf( wfs, w_cx, nevt, 232 );
    retrieveWf( wfs, w_cx, nevt, 234 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 255 );
    retrieveWf( wfs, w_cx, nevt, 355 );
    retrieveWf( wfs, w_cx, nevt, 361 );
    retrieveWf( wfs, w_cx, nevt, 362 );
    retrieveWf( wfs, w_cx, nevt, 438 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 473 );
    retrieveWf( wfs, w_cx, nevt, 479 );
    retrieveWf( wfs, w_cx, nevt, 480 );
    retrieveWf( wfs, w_cx, nevt, 481 );
    retrieveWf( wfs, w_cx, nevt, 488 );
    retrieveWf( wfs, w_cx, nevt, 510 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 519 );
    retrieveWf( wfs, w_cx, nevt, 521 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 523 );
    retrieveWf( wfs, w_cx, nevt, 524 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 530 );
    retrieveWf( wfs, w_cx, nevt, 532 );
    retrieveWf( wfs, w_cx, nevt, 533 );
    retrieveWf( wfs, w_cx, nevt, 538 );
    retrieveWf( wfs, w_cx, nevt, 539 );
    retrieveWf( wfs, w_cx, nevt, 543 );
    retrieveWf( wfs, w_cx, nevt, 547 );
    retrieveWf( wfs, w_cx, nevt, 550 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 556 );
    retrieveWf( wfs, w_cx, nevt, 557 );
    retrieveWf( wfs, w_cx, nevt, 558 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 560 );
    retrieveWf( wfs, w_cx, nevt, 562 );
    retrieveWf( wfs, w_cx, nevt, 571 );
    retrieveWf( wfs, w_cx, nevt, 577 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 580 );
    retrieveWf( wfs, w_cx, nevt, 581 );
    retrieveWf( wfs, w_cx, nevt, 585 );
    retrieveWf( wfs, w_cx, nevt, 586 );
    retrieveWf( wfs, w_cx, nevt, 588 );
    retrieveWf( wfs, w_cx, nevt, 589 );
    retrieveWf( wfs, w_cx, nevt, 590 );
    retrieveWf( wfs, w_cx, nevt, 591 );
    retrieveWf( wfs, w_cx, nevt, 592 );
    retrieveWf( wfs, w_cx, nevt, 593 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 596 );
#endif
#endif

    // *** DIAGRAM 5901 OF 15495 ***
    // Wavefunction(s) for diagram number 5901
    // (none)
    // Amplitude(s) for diagram number 5901
    FFV1_0( w_fp[202], w_fp[215], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];

    // *** DIAGRAM 5902 OF 15495 ***
    // Wavefunction(s) for diagram number 5902
    // (none)
    // Amplitude(s) for diagram number 5902
    FFV1_0( w_fp[255], w_fp[519], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5903 OF 15495 ***
    // Wavefunction(s) for diagram number 5903
    // (none)
    // Amplitude(s) for diagram number 5903
    FFV1_0( w_fp[202], w_fp[519], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5904 OF 15495 ***
    // Wavefunction(s) for diagram number 5904
    // (none)
    // Amplitude(s) for diagram number 5904
    FFV1_0( w_fp[580], w_fp[473], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5905 OF 15495 ***
    // Wavefunction(s) for diagram number 5905
    // (none)
    // Amplitude(s) for diagram number 5905
    FFV1_0( w_fp[580], w_fp[225], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5906 OF 15495 ***
    // Wavefunction(s) for diagram number 5906
    // (none)
    // Amplitude(s) for diagram number 5906
    FFV1_0( w_fp[174], w_fp[473], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];

    // *** DIAGRAM 5907 OF 15495 ***
    // Wavefunction(s) for diagram number 5907
    // (none)
    // Amplitude(s) for diagram number 5907
    FFV1_0( w_fp[255], w_fp[215], w_fp[577], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 5908 OF 15495 ***
    // Wavefunction(s) for diagram number 5908
    // (none)
    // Amplitude(s) for diagram number 5908
    VVV1_0( w_fp[577], w_fp[1], w_fp[232], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5909 OF 15495 ***
    // Wavefunction(s) for diagram number 5909
    // (none)
    // Amplitude(s) for diagram number 5909
    FFV1_0( w_fp[202], w_fp[473], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5910 OF 15495 ***
    // Wavefunction(s) for diagram number 5910
    // (none)
    // Amplitude(s) for diagram number 5910
    FFV1_0( w_fp[255], w_fp[225], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5911 OF 15495 ***
    // Wavefunction(s) for diagram number 5911
    // (none)
    // Amplitude(s) for diagram number 5911
    FFV1_0( w_fp[174], w_fp[215], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[630] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[632] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[215], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[654] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5912 OF 15495 ***
    // Wavefunction(s) for diagram number 5912
    // (none)
    // Amplitude(s) for diagram number 5912
    VVV1_0( w_fp[521], w_fp[228], w_fp[86], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 5913 OF 15495 ***
    // Wavefunction(s) for diagram number 5913
    // (none)
    // Amplitude(s) for diagram number 5913
    FFV1_0( w_fp[3], w_fp[234], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[667] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[670] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[709] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[712] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5914 OF 15495 ***
    // Wavefunction(s) for diagram number 5914
    // (none)
    // Amplitude(s) for diagram number 5914
    FFV1_0( w_fp[206], w_fp[215], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5915 OF 15495 ***
    // Wavefunction(s) for diagram number 5915
    // (none)
    // Amplitude(s) for diagram number 5915
    FFV1_0( w_fp[3], w_fp[519], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[615] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[617] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[675] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[677] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5916 OF 15495 ***
    // Wavefunction(s) for diagram number 5916
    // (none)
    // Amplitude(s) for diagram number 5916
    FFV1_0( w_fp[206], w_fp[519], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 5917 OF 15495 ***
    // Wavefunction(s) for diagram number 5917
    // (none)
    // Amplitude(s) for diagram number 5917
    FFV1_0( w_fp[530], w_fp[473], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];

    // *** DIAGRAM 5918 OF 15495 ***
    // Wavefunction(s) for diagram number 5918
    // (none)
    // Amplitude(s) for diagram number 5918
    FFV1_0( w_fp[530], w_fp[215], w_fp[362], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[634] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[668] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[669] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[710] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[711] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5919 OF 15495 ***
    // Wavefunction(s) for diagram number 5919
    // (none)
    // Amplitude(s) for diagram number 5919
    FFV1_0( w_fp[530], w_fp[234], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[668] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];

    // *** DIAGRAM 5920 OF 15495 ***
    // Wavefunction(s) for diagram number 5920
    // (none)
    // Amplitude(s) for diagram number 5920
    FFV1_0( w_fp[3], w_fp[473], w_fp[571], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[634] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[635] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[644] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[645] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5921 OF 15495 ***
    // Wavefunction(s) for diagram number 5921
    // (none)
    // Amplitude(s) for diagram number 5921
    VVV1_0( w_fp[571], w_fp[1], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] -= amp_sv[0];
    jamp_sv[617] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[634] -= amp_sv[0];
    jamp_sv[635] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[644] += amp_sv[0];
    jamp_sv[645] -= amp_sv[0];
    jamp_sv[667] += amp_sv[0];
    jamp_sv[670] -= amp_sv[0];
    jamp_sv[675] += amp_sv[0];
    jamp_sv[677] -= amp_sv[0];
    jamp_sv[709] -= amp_sv[0];
    jamp_sv[712] += amp_sv[0];

    // *** DIAGRAM 5922 OF 15495 ***
    // Wavefunction(s) for diagram number 5922
    // (none)
    // Amplitude(s) for diagram number 5922
    FFV1_0( w_fp[206], w_fp[473], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[626] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];

    // *** DIAGRAM 5923 OF 15495 ***
    // Wavefunction(s) for diagram number 5923
    // (none)
    // Amplitude(s) for diagram number 5923
    VVV1_0( w_fp[479], w_fp[362], w_fp[228], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];

    // *** DIAGRAM 5924 OF 15495 ***
    // Wavefunction(s) for diagram number 5924
    // (none)
    // Amplitude(s) for diagram number 5924
    FFV1_0( w_fp[3], w_fp[215], w_fp[594], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[668] -= amp_sv[0];
    jamp_sv[669] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[710] += amp_sv[0];
    jamp_sv[711] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[595], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[615] += amp_sv[0];
    jamp_sv[617] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[634] += amp_sv[0];
    jamp_sv[635] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[644] -= amp_sv[0];
    jamp_sv[645] += amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[675] -= amp_sv[0];
    jamp_sv[677] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[215], w_fp[596], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[612] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[667] -= amp_sv[0];
    jamp_sv[668] += amp_sv[0];
    jamp_sv[669] -= amp_sv[0];
    jamp_sv[670] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[709] += amp_sv[0];
    jamp_sv[710] -= amp_sv[0];
    jamp_sv[711] += amp_sv[0];
    jamp_sv[712] -= amp_sv[0];

    // *** DIAGRAM 5925 OF 15495 ***
    // Wavefunction(s) for diagram number 5925
    // (none)
    // Amplitude(s) for diagram number 5925
    VVVV1_0( w_fp[521], w_fp[239], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    VVVV3_0( w_fp[521], w_fp[239], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    VVVV4_0( w_fp[521], w_fp[239], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 5926 OF 15495 ***
    // Wavefunction(s) for diagram number 5926
    // (none)
    // Amplitude(s) for diagram number 5926
    VVV1_0( w_fp[239], w_fp[7], w_fp[480], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 5927 OF 15495 ***
    // Wavefunction(s) for diagram number 5927
    // (none)
    // Amplitude(s) for diagram number 5927
    VVV1_0( w_fp[239], w_fp[6], w_fp[481], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 5928 OF 15495 ***
    // Wavefunction(s) for diagram number 5928
    FFV1_1( w_fp[2], w_fp[521], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[529] );
    // Amplitude(s) for diagram number 5928
    FFV1_0( w_fp[198], w_fp[529], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];

    // *** DIAGRAM 5929 OF 15495 ***
    // Wavefunction(s) for diagram number 5929
    // (none)
    // Amplitude(s) for diagram number 5929
    FFV1_0( w_fp[198], w_fp[2], w_fp[481], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5930 OF 15495 ***
    // Wavefunction(s) for diagram number 5930
    // (none)
    // Amplitude(s) for diagram number 5930
    FFV1_0( w_fp[199], w_fp[529], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];

    // *** DIAGRAM 5931 OF 15495 ***
    // Wavefunction(s) for diagram number 5931
    // (none)
    // Amplitude(s) for diagram number 5931
    FFV1_0( w_fp[199], w_fp[2], w_fp[480], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5932 OF 15495 ***
    // Wavefunction(s) for diagram number 5932
    // (none)
    // Amplitude(s) for diagram number 5932
    FFV1_0( w_fp[355], w_fp[550], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5933 OF 15495 ***
    // Wavefunction(s) for diagram number 5933
    // (none)
    // Amplitude(s) for diagram number 5933
    FFV1_0( w_fp[355], w_fp[590], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5934 OF 15495 ***
    // Wavefunction(s) for diagram number 5934
    FFV1_1( w_fp[532], w_fp[1], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[575] );
    // Amplitude(s) for diagram number 5934
    FFV1_0( w_fp[198], w_fp[575], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5935 OF 15495 ***
    // Wavefunction(s) for diagram number 5935
    // (none)
    // Amplitude(s) for diagram number 5935
    FFV1_0( w_fp[198], w_fp[590], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5936 OF 15495 ***
    // Wavefunction(s) for diagram number 5936
    // (none)
    // Amplitude(s) for diagram number 5936
    FFV1_0( w_fp[199], w_fp[575], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5937 OF 15495 ***
    // Wavefunction(s) for diagram number 5937
    // (none)
    // Amplitude(s) for diagram number 5937
    FFV1_0( w_fp[199], w_fp[550], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5938 OF 15495 ***
    // Wavefunction(s) for diagram number 5938
    // (none)
    // Amplitude(s) for diagram number 5938
    FFV1_0( w_fp[355], w_fp[518], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];

    // *** DIAGRAM 5939 OF 15495 ***
    // Wavefunction(s) for diagram number 5939
    // (none)
    // Amplitude(s) for diagram number 5939
    FFV1_0( w_fp[355], w_fp[2], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5940 OF 15495 ***
    // Wavefunction(s) for diagram number 5940
    // (none)
    // Amplitude(s) for diagram number 5940
    VVVV1_0( w_fp[547], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    VVVV3_0( w_fp[547], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    VVVV4_0( w_fp[547], w_fp[1], w_fp[239], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 5941 OF 15495 ***
    // Wavefunction(s) for diagram number 5941
    // (none)
    // Amplitude(s) for diagram number 5941
    VVV1_0( w_fp[239], w_fp[7], w_fp[585], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    jamp_sv[643] -= amp_sv[0];
    jamp_sv[646] += amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 5942 OF 15495 ***
    // Wavefunction(s) for diagram number 5942
    // (none)
    // Amplitude(s) for diagram number 5942
    VVV1_0( w_fp[1], w_fp[239], w_fp[516], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 5943 OF 15495 ***
    // Wavefunction(s) for diagram number 5943
    // (none)
    // Amplitude(s) for diagram number 5943
    FFV1_0( w_fp[199], w_fp[2], w_fp[585], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5944 OF 15495 ***
    // Wavefunction(s) for diagram number 5944
    // (none)
    // Amplitude(s) for diagram number 5944
    FFV1_0( w_fp[199], w_fp[518], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += amp_sv[0];
    jamp_sv[373] -= amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];

    // *** DIAGRAM 5945 OF 15495 ***
    // Wavefunction(s) for diagram number 5945
    // (none)
    // Amplitude(s) for diagram number 5945
    FFV1_0( w_fp[355], w_fp[438], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];

    // *** DIAGRAM 5946 OF 15495 ***
    // Wavefunction(s) for diagram number 5946
    // (none)
    // Amplitude(s) for diagram number 5946
    FFV1_0( w_fp[355], w_fp[2], w_fp[543], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5947 OF 15495 ***
    // Wavefunction(s) for diagram number 5947
    // (none)
    // Amplitude(s) for diagram number 5947
    VVVV1_0( w_fp[488], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[1], w_fp[239], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];

    // *** DIAGRAM 5948 OF 15495 ***
    // Wavefunction(s) for diagram number 5948
    // (none)
    // Amplitude(s) for diagram number 5948
    VVV1_0( w_fp[239], w_fp[6], w_fp[562], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[523] -= amp_sv[0];
    jamp_sv[526] += amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 5949 OF 15495 ***
    // Wavefunction(s) for diagram number 5949
    // (none)
    // Amplitude(s) for diagram number 5949
    VVV1_0( w_fp[1], w_fp[239], w_fp[543], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];

    // *** DIAGRAM 5950 OF 15495 ***
    // Wavefunction(s) for diagram number 5950
    // (none)
    // Amplitude(s) for diagram number 5950
    FFV1_0( w_fp[198], w_fp[2], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5951 OF 15495 ***
    // Wavefunction(s) for diagram number 5951
    // (none)
    // Amplitude(s) for diagram number 5951
    FFV1_0( w_fp[198], w_fp[438], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += amp_sv[0];
    jamp_sv[379] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];

    // *** DIAGRAM 5952 OF 15495 ***
    // Wavefunction(s) for diagram number 5952
    // (none)
    // Amplitude(s) for diagram number 5952
    VVV1_0( w_fp[557], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[507] += amp_sv[0];
    jamp_sv[517] -= amp_sv[0];
    jamp_sv[613] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[673] += amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    VVV1_0( w_fp[556], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[195] -= amp_sv[0];
    jamp_sv[205] += amp_sv[0];
    jamp_sv[373] += amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[616] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[643] += amp_sv[0];
    jamp_sv[646] -= amp_sv[0];
    jamp_sv[676] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0( w_fp[555], w_fp[239], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 5953 OF 15495 ***
    // Wavefunction(s) for diagram number 5953
    // (none)
    // Amplitude(s) for diagram number 5953
    FFV1_0( w_fp[199], w_fp[2], w_fp[557], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[61] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[195] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[205] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[373] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[199], w_fp[2], w_fp[555], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[493] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[507] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[517] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[553] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5954 OF 15495 ***
    // Wavefunction(s) for diagram number 5954
    // (none)
    // Amplitude(s) for diagram number 5954
    VVV1_0( w_fp[560], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += amp_sv[0];
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[493] -= amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[553] += amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[627] += amp_sv[0];
    jamp_sv[637] -= amp_sv[0];
    VVV1_0( w_fp[559], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[219] -= amp_sv[0];
    jamp_sv[229] += amp_sv[0];
    jamp_sv[379] += amp_sv[0];
    jamp_sv[496] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[523] += amp_sv[0];
    jamp_sv[526] -= amp_sv[0];
    jamp_sv[556] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];
    VVV1_0( w_fp[558], w_fp[239], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[493] += amp_sv[0];
    jamp_sv[507] -= amp_sv[0];
    jamp_sv[517] += amp_sv[0];
    jamp_sv[553] -= amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[613] += amp_sv[0];
    jamp_sv[627] -= amp_sv[0];
    jamp_sv[637] += amp_sv[0];
    jamp_sv[673] -= amp_sv[0];

    // *** DIAGRAM 5955 OF 15495 ***
    // Wavefunction(s) for diagram number 5955
    // (none)
    // Amplitude(s) for diagram number 5955
    FFV1_0( w_fp[198], w_fp[2], w_fp[560], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[67] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[219] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[229] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[379] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[198], w_fp[2], w_fp[558], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[613] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[627] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[637] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[673] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5956 OF 15495 ***
    // Wavefunction(s) for diagram number 5956
    // (none)
    // Amplitude(s) for diagram number 5956
    FFV1_0( w_fp[355], w_fp[2], w_fp[538], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[527], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[355], w_fp[2], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[496] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[556] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[616] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[676] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5957 OF 15495 ***
    // Wavefunction(s) for diagram number 5957
    // (none)
    // Amplitude(s) for diagram number 5957
    VVV1_0( w_fp[538], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];
    VVV1_0( w_fp[527], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[70] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    VVV1_0( w_fp[533], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[195] += amp_sv[0];
    jamp_sv[205] -= amp_sv[0];
    jamp_sv[219] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[496] -= amp_sv[0];
    jamp_sv[556] += amp_sv[0];
    jamp_sv[616] -= amp_sv[0];
    jamp_sv[676] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 5958 OF 15495 ***
    // Wavefunction(s) for diagram number 5958
    // (none)
    // Amplitude(s) for diagram number 5958
    VVV1_0( w_fp[521], w_fp[239], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 5959 OF 15495 ***
    // Wavefunction(s) for diagram number 5959
    // (none)
    // Amplitude(s) for diagram number 5959
    FFV1_0( w_fp[196], w_fp[128], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5960 OF 15495 ***
    // Wavefunction(s) for diagram number 5960
    // (none)
    // Amplitude(s) for diagram number 5960
    FFV1_0( w_fp[129], w_fp[2], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[129] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5961 OF 15495 ***
    // Wavefunction(s) for diagram number 5961
    // (none)
    // Amplitude(s) for diagram number 5961
    FFV1_0( w_fp[355], w_fp[532], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];

    // *** DIAGRAM 5962 OF 15495 ***
    // Wavefunction(s) for diagram number 5962
    // (none)
    // Amplitude(s) for diagram number 5962
    FFV1_0( w_fp[196], w_fp[532], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[53] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[64] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[363] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[365] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5963 OF 15495 ***
    // Wavefunction(s) for diagram number 5963
    // (none)
    // Amplitude(s) for diagram number 5963
    FFV1_0( w_fp[129], w_fp[532], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];

    // *** DIAGRAM 5964 OF 15495 ***
    // Wavefunction(s) for diagram number 5964
    // (none)
    // Amplitude(s) for diagram number 5964
    FFV1_0( w_fp[453], w_fp[2], w_fp[361], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[214] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[238] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[583] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[703] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5965 OF 15495 ***
    // Wavefunction(s) for diagram number 5965
    // (none)
    // Amplitude(s) for diagram number 5965
    FFV1_0( w_fp[453], w_fp[128], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];

    // *** DIAGRAM 5966 OF 15495 ***
    // Wavefunction(s) for diagram number 5966
    // (none)
    // Amplitude(s) for diagram number 5966
    FFV1_0( w_fp[355], w_fp[2], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[70] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[376] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[382] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[594] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[714] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5967 OF 15495 ***
    // Wavefunction(s) for diagram number 5967
    // (none)
    // Amplitude(s) for diagram number 5967
    VVV1_0( w_fp[581], w_fp[1], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] -= amp_sv[0];
    jamp_sv[70] += amp_sv[0];
    jamp_sv[129] += amp_sv[0];
    jamp_sv[131] -= amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[214] += amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[238] -= amp_sv[0];
    jamp_sv[376] += amp_sv[0];
    jamp_sv[382] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 5968 OF 15495 ***
    // Wavefunction(s) for diagram number 5968
    // (none)
    // Amplitude(s) for diagram number 5968
    FFV1_0( w_fp[355], w_fp[128], w_fp[479], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[580] += amp_sv[0];
    jamp_sv[594] -= amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[714] += amp_sv[0];

    // *** DIAGRAM 5969 OF 15495 ***
    // Wavefunction(s) for diagram number 5969
    // (none)
    // Amplitude(s) for diagram number 5969
    VVV1_0( w_fp[479], w_fp[361], w_fp[239], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];

    // *** DIAGRAM 5970 OF 15495 ***
    // Wavefunction(s) for diagram number 5970
    // (none)
    // Amplitude(s) for diagram number 5970
    FFV1_0( w_fp[196], w_fp[2], w_fp[593], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[363] += amp_sv[0];
    jamp_sv[365] -= amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[583] -= amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[703] += amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[592], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[64] += amp_sv[0];
    jamp_sv[70] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[214] -= amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[238] += amp_sv[0];
    jamp_sv[376] -= amp_sv[0];
    jamp_sv[382] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];
    FFV1_0( w_fp[196], w_fp[2], w_fp[591], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[129] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[363] -= amp_sv[0];
    jamp_sv[365] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[583] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[594] += amp_sv[0];
    jamp_sv[700] += amp_sv[0];
    jamp_sv[703] -= amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[714] -= amp_sv[0];

    // *** DIAGRAM 5971 OF 15495 ***
    // Wavefunction(s) for diagram number 5971
    // (none)
    // Amplitude(s) for diagram number 5971
    VVVV1_0( w_fp[521], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    VVVV3_0( w_fp[521], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    VVVV4_0( w_fp[521], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 5972 OF 15495 ***
    // Wavefunction(s) for diagram number 5972
    // (none)
    // Amplitude(s) for diagram number 5972
    VVV1_0( w_fp[150], w_fp[7], w_fp[510], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 5973 OF 15495 ***
    // Wavefunction(s) for diagram number 5973
    // (none)
    // Amplitude(s) for diagram number 5973
    VVV1_0( w_fp[150], w_fp[4], w_fp[481], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 5974 OF 15495 ***
    // Wavefunction(s) for diagram number 5974
    // (none)
    // Amplitude(s) for diagram number 5974
    FFV1_0( w_fp[202], w_fp[529], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];

    // *** DIAGRAM 5975 OF 15495 ***
    // Wavefunction(s) for diagram number 5975
    // (none)
    // Amplitude(s) for diagram number 5975
    FFV1_0( w_fp[202], w_fp[2], w_fp[481], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[626] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5976 OF 15495 ***
    // Wavefunction(s) for diagram number 5976
    // (none)
    // Amplitude(s) for diagram number 5976
    FFV1_0( w_fp[176], w_fp[529], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[361] -= amp_sv[0];

    // *** DIAGRAM 5977 OF 15495 ***
    // Wavefunction(s) for diagram number 5977
    // (none)
    // Amplitude(s) for diagram number 5977
    FFV1_0( w_fp[176], w_fp[2], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5978 OF 15495 ***
    // Wavefunction(s) for diagram number 5978
    // (none)
    // Amplitude(s) for diagram number 5978
    FFV1_0( w_fp[255], w_fp[589], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5979 OF 15495 ***
    // Wavefunction(s) for diagram number 5979
    // (none)
    // Amplitude(s) for diagram number 5979
    FFV1_0( w_fp[255], w_fp[590], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5980 OF 15495 ***
    // Wavefunction(s) for diagram number 5980
    // (none)
    // Amplitude(s) for diagram number 5980
    FFV1_0( w_fp[202], w_fp[575], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[364] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5981 OF 15495 ***
    // Wavefunction(s) for diagram number 5981
    // (none)
    // Amplitude(s) for diagram number 5981
    FFV1_0( w_fp[202], w_fp[590], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5982 OF 15495 ***
    // Wavefunction(s) for diagram number 5982
    // (none)
    // Amplitude(s) for diagram number 5982
    FFV1_0( w_fp[176], w_fp[575], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5983 OF 15495 ***
    // Wavefunction(s) for diagram number 5983
    // (none)
    // Amplitude(s) for diagram number 5983
    FFV1_0( w_fp[176], w_fp[589], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5984 OF 15495 ***
    // Wavefunction(s) for diagram number 5984
    // (none)
    // Amplitude(s) for diagram number 5984
    FFV1_0( w_fp[255], w_fp[588], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] += amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[370] -= amp_sv[0];

    // *** DIAGRAM 5985 OF 15495 ***
    // Wavefunction(s) for diagram number 5985
    // (none)
    // Amplitude(s) for diagram number 5985
    FFV1_0( w_fp[255], w_fp[2], w_fp[539], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[370] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[660] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5986 OF 15495 ***
    // Wavefunction(s) for diagram number 5986
    // (none)
    // Amplitude(s) for diagram number 5986
    VVVV1_0( w_fp[577], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];
    VVVV3_0( w_fp[577], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    VVVV4_0( w_fp[577], w_fp[1], w_fp[150], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 5987 OF 15495 ***
    // Wavefunction(s) for diagram number 5987
    // (none)
    // Amplitude(s) for diagram number 5987
    VVV1_0( w_fp[150], w_fp[7], w_fp[586], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[145] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[630] -= amp_sv[0];
    jamp_sv[632] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[660] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 5988 OF 15495 ***
    // Wavefunction(s) for diagram number 5988
    // (none)
    // Amplitude(s) for diagram number 5988
    VVV1_0( w_fp[1], w_fp[150], w_fp[539], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[58] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[370] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 5989 OF 15495 ***
    // Wavefunction(s) for diagram number 5989
    // (none)
    // Amplitude(s) for diagram number 5989
    FFV1_0( w_fp[176], w_fp[2], w_fp[586], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5990 OF 15495 ***
    // Wavefunction(s) for diagram number 5990
    // (none)
    // Amplitude(s) for diagram number 5990
    FFV1_0( w_fp[176], w_fp[588], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[367] -= amp_sv[0];

    // *** DIAGRAM 5991 OF 15495 ***
    // Wavefunction(s) for diagram number 5991
    // (none)
    // Amplitude(s) for diagram number 5991
    FFV1_0( w_fp[255], w_fp[438], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] += amp_sv[0];
    jamp_sv[380] -= amp_sv[0];
    jamp_sv[614] -= amp_sv[0];
    jamp_sv[674] += amp_sv[0];

    // *** DIAGRAM 5992 OF 15495 ***
    // Wavefunction(s) for diagram number 5992
    // (none)
    // Amplitude(s) for diagram number 5992
    FFV1_0( w_fp[255], w_fp[2], w_fp[578], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[250] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[292] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[348] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[380] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[614] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[674] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5993 OF 15495 ***
    // Wavefunction(s) for diagram number 5993
    // (none)
    // Amplitude(s) for diagram number 5993
    VVVV1_0( w_fp[488], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    VVVV3_0( w_fp[488], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    VVVV4_0( w_fp[488], w_fp[1], w_fp[150], w_fp[4], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 5994 OF 15495 ***
    // Wavefunction(s) for diagram number 5994
    // (none)
    // Amplitude(s) for diagram number 5994
    VVV1_0( w_fp[150], w_fp[4], w_fp[562], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[228] -= amp_sv[0];
    jamp_sv[250] -= amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[282] -= amp_sv[0];
    jamp_sv[284] += amp_sv[0];
    jamp_sv[292] += amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[348] -= amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 5995 OF 15495 ***
    // Wavefunction(s) for diagram number 5995
    // (none)
    // Amplitude(s) for diagram number 5995
    VVV1_0( w_fp[1], w_fp[150], w_fp[578], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[68] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[162] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[380] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];

    // *** DIAGRAM 5996 OF 15495 ***
    // Wavefunction(s) for diagram number 5996
    // (none)
    // Amplitude(s) for diagram number 5996
    FFV1_0( w_fp[202], w_fp[2], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[378] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[672] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 5997 OF 15495 ***
    // Wavefunction(s) for diagram number 5997
    // (none)
    // Amplitude(s) for diagram number 5997
    FFV1_0( w_fp[202], w_fp[438], w_fp[1], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] += amp_sv[0];
    jamp_sv[378] -= amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[672] += amp_sv[0];

    // *** DIAGRAM 5998 OF 15495 ***
    // Wavefunction(s) for diagram number 5998
    // (none)
    // Amplitude(s) for diagram number 5998
    VVV1_0( w_fp[524], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += amp_sv[0];
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[265] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[361] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[654] -= amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[672] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    VVV1_0( w_fp[523], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[145] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[367] += amp_sv[0];
    jamp_sv[614] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[630] += amp_sv[0];
    jamp_sv[632] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[674] -= amp_sv[0];
    VVV1_0( w_fp[522], w_fp[150], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[361] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[654] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[660] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];

    // *** DIAGRAM 5999 OF 15495 ***
    // Wavefunction(s) for diagram number 5999
    // (none)
    // Amplitude(s) for diagram number 5999
    FFV1_0( w_fp[176], w_fp[2], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[523], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[55] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[367] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[176], w_fp[2], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[265] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[289] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[361] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 6000 OF 15495 ***
    // Wavefunction(s) for diagram number 6000
    // (none)
    // Amplitude(s) for diagram number 6000
    VVV1_0( w_fp[560], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] += amp_sv[0];
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[289] += amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[342] -= amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[364] -= amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[626] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    VVV1_0( w_fp[559], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[66] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[250] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[282] += amp_sv[0];
    jamp_sv[284] -= amp_sv[0];
    jamp_sv[292] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[378] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];
    VVV1_0( w_fp[558], w_fp[150], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[265] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[289] -= amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[342] += amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[348] += amp_sv[0];
    jamp_sv[364] += amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[626] -= amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[672] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 529 );
    storeWf( wfs, w_cx, nevt, 575 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
