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
  diagramgroup45( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 3 );
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 5 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 64 );
    retrieveWf( wfs, w_cx, nevt, 65 );
    retrieveWf( wfs, w_cx, nevt, 66 );
    retrieveWf( wfs, w_cx, nevt, 68 );
    retrieveWf( wfs, w_cx, nevt, 69 );
    retrieveWf( wfs, w_cx, nevt, 70 );
    retrieveWf( wfs, w_cx, nevt, 80 );
    retrieveWf( wfs, w_cx, nevt, 83 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 85 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 87 );
    retrieveWf( wfs, w_cx, nevt, 88 );
    retrieveWf( wfs, w_cx, nevt, 96 );
    retrieveWf( wfs, w_cx, nevt, 100 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 103 );
    retrieveWf( wfs, w_cx, nevt, 109 );
    retrieveWf( wfs, w_cx, nevt, 112 );
    retrieveWf( wfs, w_cx, nevt, 113 );
    retrieveWf( wfs, w_cx, nevt, 115 );
    retrieveWf( wfs, w_cx, nevt, 124 );
    retrieveWf( wfs, w_cx, nevt, 130 );
    retrieveWf( wfs, w_cx, nevt, 133 );
    retrieveWf( wfs, w_cx, nevt, 134 );
    retrieveWf( wfs, w_cx, nevt, 135 );
    retrieveWf( wfs, w_cx, nevt, 136 );
    retrieveWf( wfs, w_cx, nevt, 137 );
    retrieveWf( wfs, w_cx, nevt, 139 );
    retrieveWf( wfs, w_cx, nevt, 140 );
    retrieveWf( wfs, w_cx, nevt, 141 );
    retrieveWf( wfs, w_cx, nevt, 143 );
    retrieveWf( wfs, w_cx, nevt, 145 );
    retrieveWf( wfs, w_cx, nevt, 146 );
    retrieveWf( wfs, w_cx, nevt, 147 );
    retrieveWf( wfs, w_cx, nevt, 149 );
    retrieveWf( wfs, w_cx, nevt, 168 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 175 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 180 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 182 );
    retrieveWf( wfs, w_cx, nevt, 188 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 206 );
    retrieveWf( wfs, w_cx, nevt, 218 );
    retrieveWf( wfs, w_cx, nevt, 221 );
    retrieveWf( wfs, w_cx, nevt, 259 );
    retrieveWf( wfs, w_cx, nevt, 436 );
    retrieveWf( wfs, w_cx, nevt, 439 );
    retrieveWf( wfs, w_cx, nevt, 441 );
    retrieveWf( wfs, w_cx, nevt, 442 );
    retrieveWf( wfs, w_cx, nevt, 443 );
    retrieveWf( wfs, w_cx, nevt, 449 );
    retrieveWf( wfs, w_cx, nevt, 453 );
    retrieveWf( wfs, w_cx, nevt, 456 );
    retrieveWf( wfs, w_cx, nevt, 457 );
    retrieveWf( wfs, w_cx, nevt, 460 );
    retrieveWf( wfs, w_cx, nevt, 463 );
    retrieveWf( wfs, w_cx, nevt, 465 );
    retrieveWf( wfs, w_cx, nevt, 466 );
    retrieveWf( wfs, w_cx, nevt, 467 );
    retrieveWf( wfs, w_cx, nevt, 468 );
    retrieveWf( wfs, w_cx, nevt, 475 );
    retrieveWf( wfs, w_cx, nevt, 505 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 517 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 520 );
    retrieveWf( wfs, w_cx, nevt, 522 );
    retrieveWf( wfs, w_cx, nevt, 525 );
    retrieveWf( wfs, w_cx, nevt, 527 );
    retrieveWf( wfs, w_cx, nevt, 537 );
    retrieveWf( wfs, w_cx, nevt, 546 );
    retrieveWf( wfs, w_cx, nevt, 559 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 585 );
#endif
#endif

    // *** DIAGRAM 8801 OF 15495 ***
    // Wavefunction(s) for diagram number 8801
    // (none)
    // Amplitude(s) for diagram number 8801
    FFV1_0( w_fp[527], w_fp[467], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8802 OF 15495 ***
    // Wavefunction(s) for diagram number 8802
    // (none)
    // Amplitude(s) for diagram number 8802
    FFV1_0( w_fp[527], w_fp[259], w_fp[130], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];

    // *** DIAGRAM 8803 OF 15495 ***
    // Wavefunction(s) for diagram number 8803
    // (none)
    // Amplitude(s) for diagram number 8803
    VVV1_0( w_fp[516], w_fp[468], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8804 OF 15495 ***
    // Wavefunction(s) for diagram number 8804
    // (none)
    // Amplitude(s) for diagram number 8804
    FFV1_0( w_fp[168], w_fp[443], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += amp_sv[0];
    jamp_sv[149] -= amp_sv[0];
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[166] += amp_sv[0];

    // *** DIAGRAM 8805 OF 15495 ***
    // Wavefunction(s) for diagram number 8805
    // (none)
    // Amplitude(s) for diagram number 8805
    FFV1_0( w_fp[218], w_fp[259], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];

    // *** DIAGRAM 8806 OF 15495 ***
    // Wavefunction(s) for diagram number 8806
    // (none)
    // Amplitude(s) for diagram number 8806
    VVV1_0( w_fp[0], w_fp[468], w_fp[130], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8807 OF 15495 ***
    // Wavefunction(s) for diagram number 8807
    // (none)
    // Amplitude(s) for diagram number 8807
    FFV1_0( w_fp[112], w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[147] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8808 OF 15495 ***
    // Wavefunction(s) for diagram number 8808
    // (none)
    // Amplitude(s) for diagram number 8808
    FFV1_0( w_fp[218], w_fp[467], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8809 OF 15495 ***
    // Wavefunction(s) for diagram number 8809
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[468] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[579] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[556] );
    // Amplitude(s) for diagram number 8809
    FFV1_0( w_fp[168], w_fp[259], w_fp[468], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[579], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[556], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[147] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[149] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8810 OF 15495 ***
    // Wavefunction(s) for diagram number 8810
    // (none)
    // Amplitude(s) for diagram number 8810
    FFV1_0( w_fp[168], w_fp[449], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[449], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    FFV1_0( w_fp[168], w_fp[449], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[139] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];

    // *** DIAGRAM 8811 OF 15495 ***
    // Wavefunction(s) for diagram number 8811
    // (none)
    // Amplitude(s) for diagram number 8811
    FFV1_0( w_fp[527], w_fp[259], w_fp[145], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] += amp_sv[0];
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[236] += amp_sv[0];
    FFV1_0( w_fp[527], w_fp[259], w_fp[146], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[166] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[212] -= amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    FFV1_0( w_fp[527], w_fp[259], w_fp[147], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[160] -= amp_sv[0];
    jamp_sv[202] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 8812 OF 15495 ***
    // Wavefunction(s) for diagram number 8812
    VVV1P0_1( w_fp[0], w_fp[145], COUPs[0], 1.0, depCoup, 0., 0., w_fp[593] );
    VVV1P0_1( w_fp[0], w_fp[146], COUPs[0], 1.0, depCoup, 0., 0., w_fp[581] );
    VVV1P0_1( w_fp[0], w_fp[147], COUPs[0], 1.0, depCoup, 0., 0., w_fp[533] );
    // Amplitude(s) for diagram number 8812
    FFV1_0( w_fp[168], w_fp[259], w_fp[593], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[581], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[125] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[166] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[212] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[168], w_fp[259], w_fp[533], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[160] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[202] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[226] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[236] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8813 OF 15495 ***
    // Wavefunction(s) for diagram number 8813
    // (none)
    // Amplitude(s) for diagram number 8813
    FFV1_0( w_fp[202], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] -= amp_sv[0];

    // *** DIAGRAM 8814 OF 15495 ***
    // Wavefunction(s) for diagram number 8814
    // (none)
    // Amplitude(s) for diagram number 8814
    FFV1_0( w_fp[202], w_fp[585], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[140] -= amp_sv[0];

    // *** DIAGRAM 8815 OF 15495 ***
    // Wavefunction(s) for diagram number 8815
    // (none)
    // Amplitude(s) for diagram number 8815
    FFV1_0( w_fp[175], w_fp[546], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] -= amp_sv[0];

    // *** DIAGRAM 8816 OF 15495 ***
    // Wavefunction(s) for diagram number 8816
    // (none)
    // Amplitude(s) for diagram number 8816
    FFV1_0( w_fp[175], w_fp[585], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[138] -= amp_sv[0];

    // *** DIAGRAM 8817 OF 15495 ***
    // Wavefunction(s) for diagram number 8817
    // (none)
    // Amplitude(s) for diagram number 8817
    FFV1_0( w_fp[176], w_fp[546], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];

    // *** DIAGRAM 8818 OF 15495 ***
    // Wavefunction(s) for diagram number 8818
    // (none)
    // Amplitude(s) for diagram number 8818
    FFV1_0( w_fp[176], w_fp[442], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] -= amp_sv[0];

    // *** DIAGRAM 8819 OF 15495 ***
    // Wavefunction(s) for diagram number 8819
    FFV1_2( w_fp[174], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[594] );
    FFV1_2( w_fp[594], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[577] );
    // Amplitude(s) for diagram number 8819
    FFV1_0( w_fp[577], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= amp_sv[0];

    // *** DIAGRAM 8820 OF 15495 ***
    // Wavefunction(s) for diagram number 8820
    FFV1_2( w_fp[594], w_fp[7], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[578] );
    // Amplitude(s) for diagram number 8820
    FFV1_0( w_fp[578], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= amp_sv[0];

    // *** DIAGRAM 8821 OF 15495 ***
    // Wavefunction(s) for diagram number 8821
    FFV1_2( w_fp[594], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[595] );
    // Amplitude(s) for diagram number 8821
    FFV1_0( w_fp[595], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] -= amp_sv[0];

    // *** DIAGRAM 8822 OF 15495 ***
    // Wavefunction(s) for diagram number 8822
    // (none)
    // Amplitude(s) for diagram number 8822
    FFV1_0( w_fp[578], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] -= amp_sv[0];

    // *** DIAGRAM 8823 OF 15495 ***
    // Wavefunction(s) for diagram number 8823
    // (none)
    // Amplitude(s) for diagram number 8823
    FFV1_0( w_fp[595], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 8824 OF 15495 ***
    // Wavefunction(s) for diagram number 8824
    // (none)
    // Amplitude(s) for diagram number 8824
    FFV1_0( w_fp[577], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] -= amp_sv[0];

    // *** DIAGRAM 8825 OF 15495 ***
    // Wavefunction(s) for diagram number 8825
    // (none)
    // Amplitude(s) for diagram number 8825
    FFV1_0( w_fp[175], w_fp[520], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[148] -= amp_sv[0];

    // *** DIAGRAM 8826 OF 15495 ***
    // Wavefunction(s) for diagram number 8826
    FFV1_2( w_fp[175], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[530] );
    // Amplitude(s) for diagram number 8826
    FFV1_0( w_fp[530], w_fp[443], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] -= amp_sv[0];

    // *** DIAGRAM 8827 OF 15495 ***
    // Wavefunction(s) for diagram number 8827
    // (none)
    // Amplitude(s) for diagram number 8827
    FFV1_0( w_fp[176], w_fp[520], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] -= amp_sv[0];

    // *** DIAGRAM 8828 OF 15495 ***
    // Wavefunction(s) for diagram number 8828
    FFV1_2( w_fp[176], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[555] );
    // Amplitude(s) for diagram number 8828
    FFV1_0( w_fp[555], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] -= amp_sv[0];

    // *** DIAGRAM 8829 OF 15495 ***
    // Wavefunction(s) for diagram number 8829
    // (none)
    // Amplitude(s) for diagram number 8829
    FFV1_0( w_fp[202], w_fp[537], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[172] -= amp_sv[0];

    // *** DIAGRAM 8830 OF 15495 ***
    // Wavefunction(s) for diagram number 8830
    FFV1_2( w_fp[202], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[545] );
    // Amplitude(s) for diagram number 8830
    FFV1_0( w_fp[545], w_fp[436], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] -= amp_sv[0];

    // *** DIAGRAM 8831 OF 15495 ***
    // Wavefunction(s) for diagram number 8831
    // (none)
    // Amplitude(s) for diagram number 8831
    FFV1_0( w_fp[176], w_fp[537], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= amp_sv[0];

    // *** DIAGRAM 8832 OF 15495 ***
    // Wavefunction(s) for diagram number 8832
    // (none)
    // Amplitude(s) for diagram number 8832
    FFV1_0( w_fp[555], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[175] -= amp_sv[0];

    // *** DIAGRAM 8833 OF 15495 ***
    // Wavefunction(s) for diagram number 8833
    // (none)
    // Amplitude(s) for diagram number 8833
    FFV1_0( w_fp[202], w_fp[505], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[218] -= amp_sv[0];

    // *** DIAGRAM 8834 OF 15495 ***
    // Wavefunction(s) for diagram number 8834
    // (none)
    // Amplitude(s) for diagram number 8834
    FFV1_0( w_fp[545], w_fp[441], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[228] -= amp_sv[0];

    // *** DIAGRAM 8835 OF 15495 ***
    // Wavefunction(s) for diagram number 8835
    // (none)
    // Amplitude(s) for diagram number 8835
    FFV1_0( w_fp[175], w_fp[505], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= amp_sv[0];

    // *** DIAGRAM 8836 OF 15495 ***
    // Wavefunction(s) for diagram number 8836
    // (none)
    // Amplitude(s) for diagram number 8836
    FFV1_0( w_fp[530], w_fp[441], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[222] -= amp_sv[0];

    // *** DIAGRAM 8837 OF 15495 ***
    // Wavefunction(s) for diagram number 8837
    // (none)
    // Amplitude(s) for diagram number 8837
    FFV1_0( w_fp[149], w_fp[449], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8838 OF 15495 ***
    // Wavefunction(s) for diagram number 8838
    // (none)
    // Amplitude(s) for diagram number 8838
    FFV1_0( w_fp[176], w_fp[449], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8839 OF 15495 ***
    // Wavefunction(s) for diagram number 8839
    // (none)
    // Amplitude(s) for diagram number 8839
    FFV1_0( w_fp[174], w_fp[449], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];

    // *** DIAGRAM 8840 OF 15495 ***
    // Wavefunction(s) for diagram number 8840
    // (none)
    // Amplitude(s) for diagram number 8840
    FFV1_0( w_fp[594], w_fp[457], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8841 OF 15495 ***
    // Wavefunction(s) for diagram number 8841
    // (none)
    // Amplitude(s) for diagram number 8841
    FFV1_0( w_fp[594], w_fp[441], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8842 OF 15495 ***
    // Wavefunction(s) for diagram number 8842
    // (none)
    // Amplitude(s) for diagram number 8842
    FFV1_0( w_fp[594], w_fp[259], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 8843 OF 15495 ***
    // Wavefunction(s) for diagram number 8843
    VVV1P0_1( w_fp[0], w_fp[66], COUPs[0], 1.0, depCoup, 0., 0., w_fp[514] );
    // Amplitude(s) for diagram number 8843
    VVV1_0( w_fp[514], w_fp[475], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8844 OF 15495 ***
    // Wavefunction(s) for diagram number 8844
    // (none)
    // Amplitude(s) for diagram number 8844
    FFV1_0( w_fp[174], w_fp[441], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 8845 OF 15495 ***
    // Wavefunction(s) for diagram number 8845
    // (none)
    // Amplitude(s) for diagram number 8845
    FFV1_0( w_fp[176], w_fp[259], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];

    // *** DIAGRAM 8846 OF 15495 ***
    // Wavefunction(s) for diagram number 8846
    // (none)
    // Amplitude(s) for diagram number 8846
    VVV1_0( w_fp[0], w_fp[475], w_fp[69], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8847 OF 15495 ***
    // Wavefunction(s) for diagram number 8847
    // (none)
    // Amplitude(s) for diagram number 8847
    FFV1_0( w_fp[176], w_fp[457], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8848 OF 15495 ***
    // Wavefunction(s) for diagram number 8848
    // (none)
    // Amplitude(s) for diagram number 8848
    FFV1_0( w_fp[149], w_fp[441], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8849 OF 15495 ***
    // Wavefunction(s) for diagram number 8849
    VVVV1P0_1( w_fp[0], w_fp[66], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[521] );
    VVVV3P0_1( w_fp[0], w_fp[66], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[526] );
    VVVV4P0_1( w_fp[0], w_fp[66], w_fp[7], COUPs[2], 1.0, depCoup, 0., 0., w_fp[524] );
    // Amplitude(s) for diagram number 8849
    FFV1_0( w_fp[174], w_fp[259], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8850 OF 15495 ***
    // Wavefunction(s) for diagram number 8850
    // (none)
    // Amplitude(s) for diagram number 8850
    FFV1_0( w_fp[83], w_fp[449], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8851 OF 15495 ***
    // Wavefunction(s) for diagram number 8851
    // (none)
    // Amplitude(s) for diagram number 8851
    FFV1_0( w_fp[175], w_fp[449], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8852 OF 15495 ***
    // Wavefunction(s) for diagram number 8852
    // (none)
    // Amplitude(s) for diagram number 8852
    FFV1_0( w_fp[174], w_fp[449], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[130] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];

    // *** DIAGRAM 8853 OF 15495 ***
    // Wavefunction(s) for diagram number 8853
    // (none)
    // Amplitude(s) for diagram number 8853
    FFV1_0( w_fp[594], w_fp[463], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8854 OF 15495 ***
    // Wavefunction(s) for diagram number 8854
    // (none)
    // Amplitude(s) for diagram number 8854
    FFV1_0( w_fp[594], w_fp[436], w_fp[102], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8855 OF 15495 ***
    // Wavefunction(s) for diagram number 8855
    // (none)
    // Amplitude(s) for diagram number 8855
    FFV1_0( w_fp[594], w_fp[259], w_fp[103], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[224] -= amp_sv[0];

    // *** DIAGRAM 8856 OF 15495 ***
    // Wavefunction(s) for diagram number 8856
    // (none)
    // Amplitude(s) for diagram number 8856
    VVV1_0( w_fp[518], w_fp[475], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8857 OF 15495 ***
    // Wavefunction(s) for diagram number 8857
    // (none)
    // Amplitude(s) for diagram number 8857
    FFV1_0( w_fp[174], w_fp[436], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];

    // *** DIAGRAM 8858 OF 15495 ***
    // Wavefunction(s) for diagram number 8858
    // (none)
    // Amplitude(s) for diagram number 8858
    FFV1_0( w_fp[175], w_fp[259], w_fp[518], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[162] -= amp_sv[0];
    jamp_sv[222] += amp_sv[0];

    // *** DIAGRAM 8859 OF 15495 ***
    // Wavefunction(s) for diagram number 8859
    // (none)
    // Amplitude(s) for diagram number 8859
    VVV1_0( w_fp[0], w_fp[475], w_fp[103], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8860 OF 15495 ***
    // Wavefunction(s) for diagram number 8860
    // (none)
    // Amplitude(s) for diagram number 8860
    FFV1_0( w_fp[175], w_fp[463], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[162] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8861 OF 15495 ***
    // Wavefunction(s) for diagram number 8861
    // (none)
    // Amplitude(s) for diagram number 8861
    FFV1_0( w_fp[83], w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8862 OF 15495 ***
    // Wavefunction(s) for diagram number 8862
    VVVV1P0_1( w_fp[0], w_fp[102], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[562] );
    VVVV3P0_1( w_fp[0], w_fp[102], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[471] );
    VVVV4P0_1( w_fp[0], w_fp[102], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[582] );
    // Amplitude(s) for diagram number 8862
    FFV1_0( w_fp[174], w_fp[259], w_fp[562], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[471], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[582], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[162] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[222] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8863 OF 15495 ***
    // Wavefunction(s) for diagram number 8863
    // (none)
    // Amplitude(s) for diagram number 8863
    FFV1_0( w_fp[202], w_fp[449], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8864 OF 15495 ***
    // Wavefunction(s) for diagram number 8864
    // (none)
    // Amplitude(s) for diagram number 8864
    FFV1_0( w_fp[109], w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8865 OF 15495 ***
    // Wavefunction(s) for diagram number 8865
    // (none)
    // Amplitude(s) for diagram number 8865
    FFV1_0( w_fp[174], w_fp[449], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];

    // *** DIAGRAM 8866 OF 15495 ***
    // Wavefunction(s) for diagram number 8866
    // (none)
    // Amplitude(s) for diagram number 8866
    FFV1_0( w_fp[594], w_fp[443], w_fp[100], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8867 OF 15495 ***
    // Wavefunction(s) for diagram number 8867
    // (none)
    // Amplitude(s) for diagram number 8867
    FFV1_0( w_fp[594], w_fp[466], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8868 OF 15495 ***
    // Wavefunction(s) for diagram number 8868
    // (none)
    // Amplitude(s) for diagram number 8868
    FFV1_0( w_fp[594], w_fp[259], w_fp[124], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 8869 OF 15495 ***
    // Wavefunction(s) for diagram number 8869
    // (none)
    // Amplitude(s) for diagram number 8869
    VVV1_0( w_fp[522], w_fp[475], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8870 OF 15495 ***
    // Wavefunction(s) for diagram number 8870
    // (none)
    // Amplitude(s) for diagram number 8870
    FFV1_0( w_fp[174], w_fp[443], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += amp_sv[0];
    jamp_sv[148] -= amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[164] += amp_sv[0];

    // *** DIAGRAM 8871 OF 15495 ***
    // Wavefunction(s) for diagram number 8871
    // (none)
    // Amplitude(s) for diagram number 8871
    FFV1_0( w_fp[202], w_fp[259], w_fp[522], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];

    // *** DIAGRAM 8872 OF 15495 ***
    // Wavefunction(s) for diagram number 8872
    // (none)
    // Amplitude(s) for diagram number 8872
    VVV1_0( w_fp[0], w_fp[475], w_fp[124], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8873 OF 15495 ***
    // Wavefunction(s) for diagram number 8873
    // (none)
    // Amplitude(s) for diagram number 8873
    FFV1_0( w_fp[109], w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[145] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8874 OF 15495 ***
    // Wavefunction(s) for diagram number 8874
    // (none)
    // Amplitude(s) for diagram number 8874
    FFV1_0( w_fp[202], w_fp[466], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8875 OF 15495 ***
    // Wavefunction(s) for diagram number 8875
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[100], COUPs[2], 1.0, depCoup, 0., 0., w_fp[475] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[100], COUPs[2], 1.0, depCoup, 0., 0., w_fp[553] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[100], COUPs[2], 1.0, depCoup, 0., 0., w_fp[592] );
    // Amplitude(s) for diagram number 8875
    FFV1_0( w_fp[174], w_fp[259], w_fp[475], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[553], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[592], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[145] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[148] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[228] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8876 OF 15495 ***
    // Wavefunction(s) for diagram number 8876
    // (none)
    // Amplitude(s) for diagram number 8876
    FFV1_0( w_fp[174], w_fp[449], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[449], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    FFV1_0( w_fp[174], w_fp[449], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];

    // *** DIAGRAM 8877 OF 15495 ***
    // Wavefunction(s) for diagram number 8877
    // (none)
    // Amplitude(s) for diagram number 8877
    FFV1_0( w_fp[594], w_fp[259], w_fp[139], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] += amp_sv[0];
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    FFV1_0( w_fp[594], w_fp[259], w_fp[140], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[164] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    FFV1_0( w_fp[594], w_fp[259], w_fp[141], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 8878 OF 15495 ***
    // Wavefunction(s) for diagram number 8878
    VVV1P0_1( w_fp[0], w_fp[139], COUPs[0], 1.0, depCoup, 0., 0., w_fp[591] );
    VVV1P0_1( w_fp[0], w_fp[140], COUPs[0], 1.0, depCoup, 0., 0., w_fp[539] );
    VVV1P0_1( w_fp[0], w_fp[141], COUPs[0], 1.0, depCoup, 0., 0., w_fp[531] );
    // Amplitude(s) for diagram number 8878
    FFV1_0( w_fp[174], w_fp[259], w_fp[591], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[539], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[124] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[164] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[174], w_fp[259], w_fp[531], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[154] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8879 OF 15495 ***
    // Wavefunction(s) for diagram number 8879
    // (none)
    // Amplitude(s) for diagram number 8879
    FFV1_0( w_fp[204], w_fp[442], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] -= amp_sv[0];

    // *** DIAGRAM 8880 OF 15495 ***
    // Wavefunction(s) for diagram number 8880
    // (none)
    // Amplitude(s) for diagram number 8880
    FFV1_0( w_fp[204], w_fp[453], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[134] -= amp_sv[0];

    // *** DIAGRAM 8881 OF 15495 ***
    // Wavefunction(s) for diagram number 8881
    // (none)
    // Amplitude(s) for diagram number 8881
    FFV1_0( w_fp[180], w_fp[546], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= amp_sv[0];

    // *** DIAGRAM 8882 OF 15495 ***
    // Wavefunction(s) for diagram number 8882
    // (none)
    // Amplitude(s) for diagram number 8882
    FFV1_0( w_fp[180], w_fp[453], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];

    // *** DIAGRAM 8883 OF 15495 ***
    // Wavefunction(s) for diagram number 8883
    // (none)
    // Amplitude(s) for diagram number 8883
    FFV1_0( w_fp[181], w_fp[546], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];

    // *** DIAGRAM 8884 OF 15495 ***
    // Wavefunction(s) for diagram number 8884
    // (none)
    // Amplitude(s) for diagram number 8884
    FFV1_0( w_fp[181], w_fp[442], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] -= amp_sv[0];

    // *** DIAGRAM 8885 OF 15495 ***
    // Wavefunction(s) for diagram number 8885
    FFV1_2( w_fp[179], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[596] );
    FFV1_2( w_fp[596], w_fp[5], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[438] );
    // Amplitude(s) for diagram number 8885
    FFV1_0( w_fp[438], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= amp_sv[0];

    // *** DIAGRAM 8886 OF 15495 ***
    // Wavefunction(s) for diagram number 8886
    FFV1_2( w_fp[596], w_fp[6], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[495] );
    // Amplitude(s) for diagram number 8886
    FFV1_0( w_fp[495], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= amp_sv[0];

    // *** DIAGRAM 8887 OF 15495 ***
    // Wavefunction(s) for diagram number 8887
    FFV1_2( w_fp[596], w_fp[4], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[435] );
    // Amplitude(s) for diagram number 8887
    FFV1_0( w_fp[435], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] -= amp_sv[0];

    // *** DIAGRAM 8888 OF 15495 ***
    // Wavefunction(s) for diagram number 8888
    // (none)
    // Amplitude(s) for diagram number 8888
    FFV1_0( w_fp[495], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] -= amp_sv[0];

    // *** DIAGRAM 8889 OF 15495 ***
    // Wavefunction(s) for diagram number 8889
    // (none)
    // Amplitude(s) for diagram number 8889
    FFV1_0( w_fp[435], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[206] -= amp_sv[0];

    // *** DIAGRAM 8890 OF 15495 ***
    // Wavefunction(s) for diagram number 8890
    // (none)
    // Amplitude(s) for diagram number 8890
    FFV1_0( w_fp[438], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] -= amp_sv[0];

    // *** DIAGRAM 8891 OF 15495 ***
    // Wavefunction(s) for diagram number 8891
    // (none)
    // Amplitude(s) for diagram number 8891
    FFV1_0( w_fp[180], w_fp[520], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[146] -= amp_sv[0];

    // *** DIAGRAM 8892 OF 15495 ***
    // Wavefunction(s) for diagram number 8892
    FFV1_2( w_fp[180], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[557] );
    // Amplitude(s) for diagram number 8892
    FFV1_0( w_fp[557], w_fp[443], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] -= amp_sv[0];

    // *** DIAGRAM 8893 OF 15495 ***
    // Wavefunction(s) for diagram number 8893
    // (none)
    // Amplitude(s) for diagram number 8893
    FFV1_0( w_fp[181], w_fp[520], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] -= amp_sv[0];

    // *** DIAGRAM 8894 OF 15495 ***
    // Wavefunction(s) for diagram number 8894
    FFV1_2( w_fp[181], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[583] );
    // Amplitude(s) for diagram number 8894
    FFV1_0( w_fp[583], w_fp[443], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] -= amp_sv[0];

    // *** DIAGRAM 8895 OF 15495 ***
    // Wavefunction(s) for diagram number 8895
    // (none)
    // Amplitude(s) for diagram number 8895
    FFV1_0( w_fp[204], w_fp[537], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] -= amp_sv[0];

    // *** DIAGRAM 8896 OF 15495 ***
    // Wavefunction(s) for diagram number 8896
    FFV1_2( w_fp[204], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[584] );
    // Amplitude(s) for diagram number 8896
    FFV1_0( w_fp[584], w_fp[436], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] -= amp_sv[0];

    // *** DIAGRAM 8897 OF 15495 ***
    // Wavefunction(s) for diagram number 8897
    // (none)
    // Amplitude(s) for diagram number 8897
    FFV1_0( w_fp[181], w_fp[537], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= amp_sv[0];

    // *** DIAGRAM 8898 OF 15495 ***
    // Wavefunction(s) for diagram number 8898
    // (none)
    // Amplitude(s) for diagram number 8898
    FFV1_0( w_fp[583], w_fp[436], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] -= amp_sv[0];

    // *** DIAGRAM 8899 OF 15495 ***
    // Wavefunction(s) for diagram number 8899
    // (none)
    // Amplitude(s) for diagram number 8899
    FFV1_0( w_fp[204], w_fp[525], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[194] -= amp_sv[0];

    // *** DIAGRAM 8900 OF 15495 ***
    // Wavefunction(s) for diagram number 8900
    // (none)
    // Amplitude(s) for diagram number 8900
    FFV1_0( w_fp[584], w_fp[439], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[204] -= amp_sv[0];

    // *** DIAGRAM 8901 OF 15495 ***
    // Wavefunction(s) for diagram number 8901
    // (none)
    // Amplitude(s) for diagram number 8901
    FFV1_0( w_fp[180], w_fp[525], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= amp_sv[0];

    // *** DIAGRAM 8902 OF 15495 ***
    // Wavefunction(s) for diagram number 8902
    // (none)
    // Amplitude(s) for diagram number 8902
    FFV1_0( w_fp[557], w_fp[439], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[198] -= amp_sv[0];

    // *** DIAGRAM 8903 OF 15495 ***
    // Wavefunction(s) for diagram number 8903
    // (none)
    // Amplitude(s) for diagram number 8903
    FFV1_0( w_fp[143], w_fp[449], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8904 OF 15495 ***
    // Wavefunction(s) for diagram number 8904
    // (none)
    // Amplitude(s) for diagram number 8904
    FFV1_0( w_fp[181], w_fp[449], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8905 OF 15495 ***
    // Wavefunction(s) for diagram number 8905
    // (none)
    // Amplitude(s) for diagram number 8905
    FFV1_0( w_fp[179], w_fp[449], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];

    // *** DIAGRAM 8906 OF 15495 ***
    // Wavefunction(s) for diagram number 8906
    // (none)
    // Amplitude(s) for diagram number 8906
    FFV1_0( w_fp[596], w_fp[457], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8907 OF 15495 ***
    // Wavefunction(s) for diagram number 8907
    // (none)
    // Amplitude(s) for diagram number 8907
    FFV1_0( w_fp[596], w_fp[439], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8908 OF 15495 ***
    // Wavefunction(s) for diagram number 8908
    // (none)
    // Amplitude(s) for diagram number 8908
    FFV1_0( w_fp[596], w_fp[259], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];

    // *** DIAGRAM 8909 OF 15495 ***
    // Wavefunction(s) for diagram number 8909
    // (none)
    // Amplitude(s) for diagram number 8909
    VVV1_0( w_fp[514], w_fp[517], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8910 OF 15495 ***
    // Wavefunction(s) for diagram number 8910
    // (none)
    // Amplitude(s) for diagram number 8910
    FFV1_0( w_fp[179], w_fp[439], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];

    // *** DIAGRAM 8911 OF 15495 ***
    // Wavefunction(s) for diagram number 8911
    // (none)
    // Amplitude(s) for diagram number 8911
    FFV1_0( w_fp[181], w_fp[259], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];

    // *** DIAGRAM 8912 OF 15495 ***
    // Wavefunction(s) for diagram number 8912
    // (none)
    // Amplitude(s) for diagram number 8912
    VVV1_0( w_fp[0], w_fp[517], w_fp[68], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8913 OF 15495 ***
    // Wavefunction(s) for diagram number 8913
    // (none)
    // Amplitude(s) for diagram number 8913
    FFV1_0( w_fp[181], w_fp[457], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8914 OF 15495 ***
    // Wavefunction(s) for diagram number 8914
    // (none)
    // Amplitude(s) for diagram number 8914
    FFV1_0( w_fp[143], w_fp[439], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8915 OF 15495 ***
    // Wavefunction(s) for diagram number 8915
    VVVV1P0_1( w_fp[0], w_fp[66], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[519] );
    VVVV3P0_1( w_fp[0], w_fp[66], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[454] );
    VVVV4P0_1( w_fp[0], w_fp[66], w_fp[6], COUPs[2], 1.0, depCoup, 0., 0., w_fp[575] );
    // Amplitude(s) for diagram number 8915
    FFV1_0( w_fp[179], w_fp[259], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8916 OF 15495 ***
    // Wavefunction(s) for diagram number 8916
    // (none)
    // Amplitude(s) for diagram number 8916
    FFV1_0( w_fp[80], w_fp[449], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8917 OF 15495 ***
    // Wavefunction(s) for diagram number 8917
    // (none)
    // Amplitude(s) for diagram number 8917
    FFV1_0( w_fp[180], w_fp[449], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8918 OF 15495 ***
    // Wavefunction(s) for diagram number 8918
    // (none)
    // Amplitude(s) for diagram number 8918
    FFV1_0( w_fp[179], w_fp[449], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[128] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];

    // *** DIAGRAM 8919 OF 15495 ***
    // Wavefunction(s) for diagram number 8919
    // (none)
    // Amplitude(s) for diagram number 8919
    FFV1_0( w_fp[596], w_fp[460], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8920 OF 15495 ***
    // Wavefunction(s) for diagram number 8920
    // (none)
    // Amplitude(s) for diagram number 8920
    FFV1_0( w_fp[596], w_fp[436], w_fp[86], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8921 OF 15495 ***
    // Wavefunction(s) for diagram number 8921
    // (none)
    // Amplitude(s) for diagram number 8921
    FFV1_0( w_fp[596], w_fp[259], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[200] -= amp_sv[0];

    // *** DIAGRAM 8922 OF 15495 ***
    // Wavefunction(s) for diagram number 8922
    // (none)
    // Amplitude(s) for diagram number 8922
    VVV1_0( w_fp[576], w_fp[517], w_fp[5], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8923 OF 15495 ***
    // Wavefunction(s) for diagram number 8923
    // (none)
    // Amplitude(s) for diagram number 8923
    FFV1_0( w_fp[179], w_fp[436], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];

    // *** DIAGRAM 8924 OF 15495 ***
    // Wavefunction(s) for diagram number 8924
    // (none)
    // Amplitude(s) for diagram number 8924
    FFV1_0( w_fp[180], w_fp[259], w_fp[576], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];

    // *** DIAGRAM 8925 OF 15495 ***
    // Wavefunction(s) for diagram number 8925
    // (none)
    // Amplitude(s) for diagram number 8925
    VVV1_0( w_fp[0], w_fp[517], w_fp[87], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8926 OF 15495 ***
    // Wavefunction(s) for diagram number 8926
    // (none)
    // Amplitude(s) for diagram number 8926
    FFV1_0( w_fp[180], w_fp[460], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[156] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8927 OF 15495 ***
    // Wavefunction(s) for diagram number 8927
    // (none)
    // Amplitude(s) for diagram number 8927
    FFV1_0( w_fp[80], w_fp[436], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8928 OF 15495 ***
    // Wavefunction(s) for diagram number 8928
    VVVV1P0_1( w_fp[0], w_fp[86], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[561] );
    VVVV3P0_1( w_fp[0], w_fp[86], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[447] );
    VVVV4P0_1( w_fp[0], w_fp[86], w_fp[5], COUPs[2], 1.0, depCoup, 0., 0., w_fp[590] );
    // Amplitude(s) for diagram number 8928
    FFV1_0( w_fp[179], w_fp[259], w_fp[561], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[447], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[590], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[156] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[198] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8929 OF 15495 ***
    // Wavefunction(s) for diagram number 8929
    // (none)
    // Amplitude(s) for diagram number 8929
    FFV1_0( w_fp[204], w_fp[449], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8930 OF 15495 ***
    // Wavefunction(s) for diagram number 8930
    // (none)
    // Amplitude(s) for diagram number 8930
    FFV1_0( w_fp[96], w_fp[449], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8931 OF 15495 ***
    // Wavefunction(s) for diagram number 8931
    // (none)
    // Amplitude(s) for diagram number 8931
    FFV1_0( w_fp[179], w_fp[449], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];

    // *** DIAGRAM 8932 OF 15495 ***
    // Wavefunction(s) for diagram number 8932
    // (none)
    // Amplitude(s) for diagram number 8932
    FFV1_0( w_fp[596], w_fp[443], w_fp[113], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8933 OF 15495 ***
    // Wavefunction(s) for diagram number 8933
    // (none)
    // Amplitude(s) for diagram number 8933
    FFV1_0( w_fp[596], w_fp[465], w_fp[4], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8934 OF 15495 ***
    // Wavefunction(s) for diagram number 8934
    // (none)
    // Amplitude(s) for diagram number 8934
    FFV1_0( w_fp[596], w_fp[259], w_fp[115], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];

    // *** DIAGRAM 8935 OF 15495 ***
    // Wavefunction(s) for diagram number 8935
    // (none)
    // Amplitude(s) for diagram number 8935
    VVV1_0( w_fp[559], w_fp[517], w_fp[4], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8936 OF 15495 ***
    // Wavefunction(s) for diagram number 8936
    // (none)
    // Amplitude(s) for diagram number 8936
    FFV1_0( w_fp[179], w_fp[443], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += amp_sv[0];
    jamp_sv[146] -= amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[158] += amp_sv[0];

    // *** DIAGRAM 8937 OF 15495 ***
    // Wavefunction(s) for diagram number 8937
    // (none)
    // Amplitude(s) for diagram number 8937
    FFV1_0( w_fp[204], w_fp[259], w_fp[559], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[204] += amp_sv[0];

    // *** DIAGRAM 8938 OF 15495 ***
    // Wavefunction(s) for diagram number 8938
    // (none)
    // Amplitude(s) for diagram number 8938
    VVV1_0( w_fp[0], w_fp[517], w_fp[115], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8939 OF 15495 ***
    // Wavefunction(s) for diagram number 8939
    // (none)
    // Amplitude(s) for diagram number 8939
    FFV1_0( w_fp[96], w_fp[443], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[144] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8940 OF 15495 ***
    // Wavefunction(s) for diagram number 8940
    // (none)
    // Amplitude(s) for diagram number 8940
    FFV1_0( w_fp[204], w_fp[465], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8941 OF 15495 ***
    // Wavefunction(s) for diagram number 8941
    VVVV1P0_1( w_fp[0], w_fp[4], w_fp[113], COUPs[2], 1.0, depCoup, 0., 0., w_fp[517] );
    VVVV3P0_1( w_fp[0], w_fp[4], w_fp[113], COUPs[2], 1.0, depCoup, 0., 0., w_fp[472] );
    VVVV4P0_1( w_fp[0], w_fp[4], w_fp[113], COUPs[2], 1.0, depCoup, 0., 0., w_fp[445] );
    // Amplitude(s) for diagram number 8941
    FFV1_0( w_fp[179], w_fp[259], w_fp[517], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[472], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[445], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[144] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[146] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[204] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8942 OF 15495 ***
    // Wavefunction(s) for diagram number 8942
    // (none)
    // Amplitude(s) for diagram number 8942
    FFV1_0( w_fp[179], w_fp[449], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[449], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    FFV1_0( w_fp[179], w_fp[449], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];

    // *** DIAGRAM 8943 OF 15495 ***
    // Wavefunction(s) for diagram number 8943
    // (none)
    // Amplitude(s) for diagram number 8943
    FFV1_0( w_fp[596], w_fp[259], w_fp[133], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] += amp_sv[0];
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    FFV1_0( w_fp[596], w_fp[259], w_fp[134], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[158] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    FFV1_0( w_fp[596], w_fp[259], w_fp[135], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];

    // *** DIAGRAM 8944 OF 15495 ***
    // Wavefunction(s) for diagram number 8944
    VVV1P0_1( w_fp[0], w_fp[133], COUPs[0], 1.0, depCoup, 0., 0., w_fp[515] );
    VVV1P0_1( w_fp[0], w_fp[134], COUPs[0], 1.0, depCoup, 0., 0., w_fp[589] );
    VVV1P0_1( w_fp[0], w_fp[135], COUPs[0], 1.0, depCoup, 0., 0., w_fp[510] );
    // Amplitude(s) for diagram number 8944
    FFV1_0( w_fp[179], w_fp[259], w_fp[515], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[589], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[158] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[179], w_fp[259], w_fp[510], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[152] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8945 OF 15495 ***
    // Wavefunction(s) for diagram number 8945
    // (none)
    // Amplitude(s) for diagram number 8945
    FFV1_0( w_fp[221], w_fp[453], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8946 OF 15495 ***
    // Wavefunction(s) for diagram number 8946
    // (none)
    // Amplitude(s) for diagram number 8946
    FFV1_0( w_fp[221], w_fp[585], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8947 OF 15495 ***
    // Wavefunction(s) for diagram number 8947
    FFV1P0_3( w_fp[3], w_fp[449], COUPs[1], 1.0, depCoup, 0., 0., w_fp[587] );
    // Amplitude(s) for diagram number 8947
    VVV1_0( w_fp[68], w_fp[7], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8948 OF 15495 ***
    // Wavefunction(s) for diagram number 8948
    // (none)
    // Amplitude(s) for diagram number 8948
    FFV1_0( w_fp[3], w_fp[585], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 8949 OF 15495 ***
    // Wavefunction(s) for diagram number 8949
    // (none)
    // Amplitude(s) for diagram number 8949
    VVV1_0( w_fp[69], w_fp[6], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8950 OF 15495 ***
    // Wavefunction(s) for diagram number 8950
    // (none)
    // Amplitude(s) for diagram number 8950
    FFV1_0( w_fp[3], w_fp[453], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];

    // *** DIAGRAM 8951 OF 15495 ***
    // Wavefunction(s) for diagram number 8951
    // (none)
    // Amplitude(s) for diagram number 8951
    FFV1_0( w_fp[3], w_fp[449], w_fp[137], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[449], w_fp[136], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[449], w_fp[182], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[134] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[138] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8952 OF 15495 ***
    // Wavefunction(s) for diagram number 8952
    // (none)
    // Amplitude(s) for diagram number 8952
    VVVV1_0( w_fp[514], w_fp[456], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVVV3_0( w_fp[514], w_fp[456], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVVV4_0( w_fp[514], w_fp[456], w_fp[6], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 8953 OF 15495 ***
    // Wavefunction(s) for diagram number 8953
    VVV1P0_1( w_fp[514], w_fp[6], COUPs[0], 1.0, depCoup, 0., 0., w_fp[597] );
    // Amplitude(s) for diagram number 8953
    VVV1_0( w_fp[456], w_fp[7], w_fp[597], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8954 OF 15495 ***
    // Wavefunction(s) for diagram number 8954
    VVV1P0_1( w_fp[514], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[598] );
    // Amplitude(s) for diagram number 8954
    VVV1_0( w_fp[456], w_fp[6], w_fp[598], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 8955 OF 15495 ***
    // Wavefunction(s) for diagram number 8955
    FFV1_2( w_fp[3], w_fp[514], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[599] );
    // Amplitude(s) for diagram number 8955
    FFV1_0( w_fp[599], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];

    // *** DIAGRAM 8956 OF 15495 ***
    // Wavefunction(s) for diagram number 8956
    // (none)
    // Amplitude(s) for diagram number 8956
    FFV1_0( w_fp[3], w_fp[439], w_fp[598], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8957 OF 15495 ***
    // Wavefunction(s) for diagram number 8957
    // (none)
    // Amplitude(s) for diagram number 8957
    FFV1_0( w_fp[599], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8958 OF 15495 ***
    // Wavefunction(s) for diagram number 8958
    // (none)
    // Amplitude(s) for diagram number 8958
    FFV1_0( w_fp[3], w_fp[441], w_fp[597], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8959 OF 15495 ***
    // Wavefunction(s) for diagram number 8959
    // (none)
    // Amplitude(s) for diagram number 8959
    VVVV1_0( w_fp[0], w_fp[456], w_fp[68], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[456], w_fp[68], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[456], w_fp[68], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];

    // *** DIAGRAM 8960 OF 15495 ***
    // Wavefunction(s) for diagram number 8960
    VVV1P0_1( w_fp[0], w_fp[456], COUPs[0], 1.0, depCoup, 0., 0., w_fp[600] );
    // Amplitude(s) for diagram number 8960
    VVV1_0( w_fp[68], w_fp[7], w_fp[600], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8961 OF 15495 ***
    // Wavefunction(s) for diagram number 8961
    VVV1P0_1( w_fp[0], w_fp[68], COUPs[0], 1.0, depCoup, 0., 0., w_fp[601] );
    // Amplitude(s) for diagram number 8961
    VVV1_0( w_fp[456], w_fp[7], w_fp[601], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8962 OF 15495 ***
    // Wavefunction(s) for diagram number 8962
    // (none)
    // Amplitude(s) for diagram number 8962
    VVVV1_0( w_fp[0], w_fp[456], w_fp[69], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    VVVV3_0( w_fp[0], w_fp[456], w_fp[69], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    VVVV4_0( w_fp[0], w_fp[456], w_fp[69], w_fp[6], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];

    // *** DIAGRAM 8963 OF 15495 ***
    // Wavefunction(s) for diagram number 8963
    // (none)
    // Amplitude(s) for diagram number 8963
    VVV1_0( w_fp[69], w_fp[6], w_fp[600], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];

    // *** DIAGRAM 8964 OF 15495 ***
    // Wavefunction(s) for diagram number 8964
    VVV1P0_1( w_fp[0], w_fp[69], COUPs[0], 1.0, depCoup, 0., 0., w_fp[602] );
    // Amplitude(s) for diagram number 8964
    VVV1_0( w_fp[456], w_fp[6], w_fp[602], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];

    // *** DIAGRAM 8965 OF 15495 ***
    // Wavefunction(s) for diagram number 8965
    // (none)
    // Amplitude(s) for diagram number 8965
    VVV1_0( w_fp[0], w_fp[456], w_fp[137], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[456], w_fp[136], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    VVV1_0( w_fp[0], w_fp[456], w_fp[182], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[132] += amp_sv[0];
    jamp_sv[134] -= amp_sv[0];
    jamp_sv[138] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[201] -= amp_sv[0];
    jamp_sv[207] += amp_sv[0];
    jamp_sv[225] -= amp_sv[0];
    jamp_sv[231] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 8966 OF 15495 ***
    // Wavefunction(s) for diagram number 8966
    // (none)
    // Amplitude(s) for diagram number 8966
    FFV1_0( w_fp[221], w_fp[525], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8967 OF 15495 ***
    // Wavefunction(s) for diagram number 8967
    FFV1_2( w_fp[221], w_fp[0], COUPs[1], 1.0, depCoup, cIPD[0], cIPD[1], w_fp[603] );
    // Amplitude(s) for diagram number 8967
    FFV1_0( w_fp[603], w_fp[439], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8968 OF 15495 ***
    // Wavefunction(s) for diagram number 8968
    // (none)
    // Amplitude(s) for diagram number 8968
    FFV1_0( w_fp[3], w_fp[525], w_fp[69], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[196] -= amp_sv[0];
    jamp_sv[197] += amp_sv[0];

    // *** DIAGRAM 8969 OF 15495 ***
    // Wavefunction(s) for diagram number 8969
    // (none)
    // Amplitude(s) for diagram number 8969
    FFV1_0( w_fp[3], w_fp[439], w_fp[602], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8970 OF 15495 ***
    // Wavefunction(s) for diagram number 8970
    // (none)
    // Amplitude(s) for diagram number 8970
    FFV1_0( w_fp[221], w_fp[505], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8971 OF 15495 ***
    // Wavefunction(s) for diagram number 8971
    // (none)
    // Amplitude(s) for diagram number 8971
    FFV1_0( w_fp[603], w_fp[441], w_fp[6], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8972 OF 15495 ***
    // Wavefunction(s) for diagram number 8972
    // (none)
    // Amplitude(s) for diagram number 8972
    FFV1_0( w_fp[3], w_fp[505], w_fp[68], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[220] -= amp_sv[0];
    jamp_sv[221] += amp_sv[0];

    // *** DIAGRAM 8973 OF 15495 ***
    // Wavefunction(s) for diagram number 8973
    // (none)
    // Amplitude(s) for diagram number 8973
    FFV1_0( w_fp[3], w_fp[441], w_fp[601], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8974 OF 15495 ***
    // Wavefunction(s) for diagram number 8974
    // (none)
    // Amplitude(s) for diagram number 8974
    VVV1_0( w_fp[519], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[200] += amp_sv[0];
    jamp_sv[206] -= amp_sv[0];
    jamp_sv[216] -= amp_sv[0];
    jamp_sv[218] += amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];
    VVV1_0( w_fp[454], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[134] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[152] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[221] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[225] += amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    VVV1_0( w_fp[575], w_fp[456], w_fp[7], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[150] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];

    // *** DIAGRAM 8975 OF 15495 ***
    // Wavefunction(s) for diagram number 8975
    // (none)
    // Amplitude(s) for diagram number 8975
    FFV1_0( w_fp[3], w_fp[441], w_fp[519], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[454], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[220] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[221] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[225] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[231] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[441], w_fp[575], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[216] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[218] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[224] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[230] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8976 OF 15495 ***
    // Wavefunction(s) for diagram number 8976
    // (none)
    // Amplitude(s) for diagram number 8976
    VVV1_0( w_fp[521], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[192] -= amp_sv[0];
    jamp_sv[194] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[224] += amp_sv[0];
    jamp_sv[230] -= amp_sv[0];
    VVV1_0( w_fp[526], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[138] -= amp_sv[0];
    jamp_sv[140] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[154] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[196] += amp_sv[0];
    jamp_sv[197] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[201] += amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[207] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    VVV1_0( w_fp[524], w_fp[456], w_fp[6], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[192] += amp_sv[0];
    jamp_sv[194] -= amp_sv[0];
    jamp_sv[200] -= amp_sv[0];
    jamp_sv[206] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[216] += amp_sv[0];
    jamp_sv[218] -= amp_sv[0];
    jamp_sv[224] -= amp_sv[0];
    jamp_sv[230] += amp_sv[0];

    // *** DIAGRAM 8977 OF 15495 ***
    // Wavefunction(s) for diagram number 8977
    // (none)
    // Amplitude(s) for diagram number 8977
    FFV1_0( w_fp[3], w_fp[439], w_fp[521], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[439], w_fp[526], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[196] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[197] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[201] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[207] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[439], w_fp[524], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[192] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[194] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[200] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[206] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8978 OF 15495 ***
    // Wavefunction(s) for diagram number 8978
    // (none)
    // Amplitude(s) for diagram number 8978
    FFV1_0( w_fp[221], w_fp[449], w_fp[84], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];

    // *** DIAGRAM 8979 OF 15495 ***
    // Wavefunction(s) for diagram number 8979
    // (none)
    // Amplitude(s) for diagram number 8979
    FFV1_0( w_fp[188], w_fp[449], w_fp[66], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];

    // *** DIAGRAM 8980 OF 15495 ***
    // Wavefunction(s) for diagram number 8980
    // (none)
    // Amplitude(s) for diagram number 8980
    FFV1_0( w_fp[3], w_fp[449], w_fp[65], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[136] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8981 OF 15495 ***
    // Wavefunction(s) for diagram number 8981
    // (none)
    // Amplitude(s) for diagram number 8981
    VVV1_0( w_fp[514], w_fp[456], w_fp[84], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8982 OF 15495 ***
    // Wavefunction(s) for diagram number 8982
    // (none)
    // Amplitude(s) for diagram number 8982
    FFV1_0( w_fp[3], w_fp[467], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[213] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[215] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[237] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[239] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8983 OF 15495 ***
    // Wavefunction(s) for diagram number 8983
    // (none)
    // Amplitude(s) for diagram number 8983
    FFV1_0( w_fp[188], w_fp[259], w_fp[514], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[121] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[127] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[150] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8984 OF 15495 ***
    // Wavefunction(s) for diagram number 8984
    // (none)
    // Amplitude(s) for diagram number 8984
    VVV1_0( w_fp[516], w_fp[456], w_fp[66], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];

    // *** DIAGRAM 8985 OF 15495 ***
    // Wavefunction(s) for diagram number 8985
    // (none)
    // Amplitude(s) for diagram number 8985
    FFV1_0( w_fp[3], w_fp[457], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[151] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[153] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[155] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[177] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[179] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8986 OF 15495 ***
    // Wavefunction(s) for diagram number 8986
    // (none)
    // Amplitude(s) for diagram number 8986
    FFV1_0( w_fp[221], w_fp[259], w_fp[516], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[137] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[143] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[210] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[211] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[234] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[235] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8987 OF 15495 ***
    // Wavefunction(s) for diagram number 8987
    // (none)
    // Amplitude(s) for diagram number 8987
    VVV1_0( w_fp[0], w_fp[456], w_fp[65], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[136] -= amp_sv[0];
    jamp_sv[137] += amp_sv[0];
    jamp_sv[142] += amp_sv[0];
    jamp_sv[143] -= amp_sv[0];
    jamp_sv[153] -= amp_sv[0];
    jamp_sv[155] += amp_sv[0];
    jamp_sv[177] += amp_sv[0];
    jamp_sv[179] -= amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8988 OF 15495 ***
    // Wavefunction(s) for diagram number 8988
    // (none)
    // Amplitude(s) for diagram number 8988
    FFV1_0( w_fp[188], w_fp[457], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[150] += amp_sv[0];
    jamp_sv[151] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[175] += amp_sv[0];

    // *** DIAGRAM 8989 OF 15495 ***
    // Wavefunction(s) for diagram number 8989
    // (none)
    // Amplitude(s) for diagram number 8989
    FFV1_0( w_fp[221], w_fp[467], w_fp[0], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[210] += amp_sv[0];
    jamp_sv[211] -= amp_sv[0];
    jamp_sv[234] -= amp_sv[0];
    jamp_sv[235] += amp_sv[0];

    // *** DIAGRAM 8990 OF 15495 ***
    // Wavefunction(s) for diagram number 8990
    VVVV1P0_1( w_fp[0], w_fp[66], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[467] );
    VVVV3P0_1( w_fp[0], w_fp[66], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[457] );
    VVVV4P0_1( w_fp[0], w_fp[66], w_fp[84], COUPs[2], 1.0, depCoup, 0., 0., w_fp[604] );
    // Amplitude(s) for diagram number 8990
    FFV1_0( w_fp[3], w_fp[259], w_fp[467], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] -= amp_sv[0];
    jamp_sv[121] += amp_sv[0];
    jamp_sv[126] += amp_sv[0];
    jamp_sv[127] -= amp_sv[0];
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[213] -= amp_sv[0];
    jamp_sv[215] += amp_sv[0];
    jamp_sv[237] += amp_sv[0];
    jamp_sv[239] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[457], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[136] += amp_sv[0];
    jamp_sv[137] -= amp_sv[0];
    jamp_sv[142] -= amp_sv[0];
    jamp_sv[143] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[153] += amp_sv[0];
    jamp_sv[155] -= amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[177] -= amp_sv[0];
    jamp_sv[179] += amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    FFV1_0( w_fp[3], w_fp[259], w_fp[604], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[120] += amp_sv[0];
    jamp_sv[121] -= amp_sv[0];
    jamp_sv[126] -= amp_sv[0];
    jamp_sv[127] += amp_sv[0];
    jamp_sv[150] -= amp_sv[0];
    jamp_sv[151] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[210] -= amp_sv[0];
    jamp_sv[211] += amp_sv[0];
    jamp_sv[213] += amp_sv[0];
    jamp_sv[215] -= amp_sv[0];
    jamp_sv[234] += amp_sv[0];
    jamp_sv[235] -= amp_sv[0];
    jamp_sv[237] -= amp_sv[0];
    jamp_sv[239] += amp_sv[0];

    // *** DIAGRAM 8991 OF 15495 ***
    // Wavefunction(s) for diagram number 8991
    // (none)
    // Amplitude(s) for diagram number 8991
    FFV1_0( w_fp[206], w_fp[442], w_fp[7], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8992 OF 15495 ***
    // Wavefunction(s) for diagram number 8992
    // (none)
    // Amplitude(s) for diagram number 8992
    FFV1_0( w_fp[206], w_fp[585], w_fp[5], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8993 OF 15495 ***
    // Wavefunction(s) for diagram number 8993
    // (none)
    // Amplitude(s) for diagram number 8993
    VVV1_0( w_fp[87], w_fp[7], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8994 OF 15495 ***
    // Wavefunction(s) for diagram number 8994
    // (none)
    // Amplitude(s) for diagram number 8994
    FFV1_0( w_fp[3], w_fp[585], w_fp[87], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[139] += amp_sv[0];
    jamp_sv[140] -= amp_sv[0];
    jamp_sv[141] += amp_sv[0];
    jamp_sv[142] -= amp_sv[0];

    // *** DIAGRAM 8995 OF 15495 ***
    // Wavefunction(s) for diagram number 8995
    // (none)
    // Amplitude(s) for diagram number 8995
    VVV1_0( w_fp[88], w_fp[5], w_fp[587], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8996 OF 15495 ***
    // Wavefunction(s) for diagram number 8996
    // (none)
    // Amplitude(s) for diagram number 8996
    FFV1_0( w_fp[3], w_fp[442], w_fp[88], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[126] += amp_sv[0];
    jamp_sv[128] -= amp_sv[0];
    jamp_sv[130] -= amp_sv[0];
    jamp_sv[131] += amp_sv[0];

    // *** DIAGRAM 8997 OF 15495 ***
    // Wavefunction(s) for diagram number 8997
    // (none)
    // Amplitude(s) for diagram number 8997
    FFV1_0( w_fp[3], w_fp[449], w_fp[85], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[449], w_fp[70], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[130] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[131] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[133] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0( w_fp[3], w_fp[449], w_fp[64], COUPs[1], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[126] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[128] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[132] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[139] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[140] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[141] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[142] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 8998 OF 15495 ***
    // Wavefunction(s) for diagram number 8998
    // (none)
    // Amplitude(s) for diagram number 8998
    VVVV1_0( w_fp[576], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[123] -= amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[133] += amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[157] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[187] += amp_sv[0];
    jamp_sv[189] += amp_sv[0];
    jamp_sv[191] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[199] -= amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    VVVV3_0( w_fp[576], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];
    VVVV4_0( w_fp[576], w_fp[456], w_fp[5], w_fp[7], COUPs[2], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 8999 OF 15495 ***
    // Wavefunction(s) for diagram number 8999
    VVV1P0_1( w_fp[576], w_fp[5], COUPs[0], 1.0, depCoup, 0., 0., w_fp[605] );
    // Amplitude(s) for diagram number 8999
    VVV1_0( w_fp[456], w_fp[7], w_fp[605], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[122] += amp_sv[0];
    jamp_sv[132] -= amp_sv[0];
    jamp_sv[156] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[198] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[228] += amp_sv[0];
    jamp_sv[229] -= amp_sv[0];
    jamp_sv[231] -= amp_sv[0];
    jamp_sv[233] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

    // *** DIAGRAM 9000 OF 15495 ***
    // Wavefunction(s) for diagram number 9000
    VVV1P0_1( w_fp[576], w_fp[7], COUPs[0], 1.0, depCoup, 0., 0., w_fp[606] );
    // Amplitude(s) for diagram number 9000
    VVV1_0( w_fp[456], w_fp[5], w_fp[606], COUPs[0], 1.0, depCoup, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[123] += amp_sv[0];
    jamp_sv[133] -= amp_sv[0];
    jamp_sv[157] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[187] -= amp_sv[0];
    jamp_sv[189] -= amp_sv[0];
    jamp_sv[191] += amp_sv[0];
    jamp_sv[199] += amp_sv[0];
    jamp_sv[217] -= amp_sv[0];
    jamp_sv[220] += amp_sv[0];
    jamp_sv[226] += amp_sv[0];
    jamp_sv[236] -= amp_sv[0];

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
    storeWf( wfs, w_cx, nevt, 435 );
    storeWf( wfs, w_cx, nevt, 438 );
    storeWf( wfs, w_cx, nevt, 445 );
    storeWf( wfs, w_cx, nevt, 447 );
    storeWf( wfs, w_cx, nevt, 454 );
    storeWf( wfs, w_cx, nevt, 457 );
    storeWf( wfs, w_cx, nevt, 467 );
    storeWf( wfs, w_cx, nevt, 468 );
    storeWf( wfs, w_cx, nevt, 471 );
    storeWf( wfs, w_cx, nevt, 472 );
    storeWf( wfs, w_cx, nevt, 475 );
    storeWf( wfs, w_cx, nevt, 495 );
    storeWf( wfs, w_cx, nevt, 510 );
    storeWf( wfs, w_cx, nevt, 514 );
    storeWf( wfs, w_cx, nevt, 515 );
    storeWf( wfs, w_cx, nevt, 517 );
    storeWf( wfs, w_cx, nevt, 519 );
    storeWf( wfs, w_cx, nevt, 521 );
    storeWf( wfs, w_cx, nevt, 524 );
    storeWf( wfs, w_cx, nevt, 526 );
    storeWf( wfs, w_cx, nevt, 530 );
    storeWf( wfs, w_cx, nevt, 531 );
    storeWf( wfs, w_cx, nevt, 533 );
    storeWf( wfs, w_cx, nevt, 539 );
    storeWf( wfs, w_cx, nevt, 545 );
    storeWf( wfs, w_cx, nevt, 553 );
    storeWf( wfs, w_cx, nevt, 555 );
    storeWf( wfs, w_cx, nevt, 556 );
    storeWf( wfs, w_cx, nevt, 557 );
    storeWf( wfs, w_cx, nevt, 561 );
    storeWf( wfs, w_cx, nevt, 562 );
    storeWf( wfs, w_cx, nevt, 575 );
    storeWf( wfs, w_cx, nevt, 577 );
    storeWf( wfs, w_cx, nevt, 578 );
    storeWf( wfs, w_cx, nevt, 579 );
    storeWf( wfs, w_cx, nevt, 581 );
    storeWf( wfs, w_cx, nevt, 582 );
    storeWf( wfs, w_cx, nevt, 583 );
    storeWf( wfs, w_cx, nevt, 584 );
    storeWf( wfs, w_cx, nevt, 587 );
    storeWf( wfs, w_cx, nevt, 589 );
    storeWf( wfs, w_cx, nevt, 590 );
    storeWf( wfs, w_cx, nevt, 591 );
    storeWf( wfs, w_cx, nevt, 592 );
    storeWf( wfs, w_cx, nevt, 593 );
    storeWf( wfs, w_cx, nevt, 594 );
    storeWf( wfs, w_cx, nevt, 595 );
    storeWf( wfs, w_cx, nevt, 596 );
    storeWf( wfs, w_cx, nevt, 597 );
    storeWf( wfs, w_cx, nevt, 598 );
    storeWf( wfs, w_cx, nevt, 599 );
    storeWf( wfs, w_cx, nevt, 600 );
    storeWf( wfs, w_cx, nevt, 601 );
    storeWf( wfs, w_cx, nevt, 602 );
    storeWf( wfs, w_cx, nevt, 603 );
    storeWf( wfs, w_cx, nevt, 604 );
    storeWf( wfs, w_cx, nevt, 605 );
    storeWf( wfs, w_cx, nevt, 606 );
#endif
#endif
  }

  //--------------------------------------------------------------------------
}
