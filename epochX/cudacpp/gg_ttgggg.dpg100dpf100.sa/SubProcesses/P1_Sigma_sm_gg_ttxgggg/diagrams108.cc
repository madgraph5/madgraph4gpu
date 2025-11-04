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
  diagramgroup108( fptype* wfs,                    // input/output wavefunctions[nwf*2*nw6*nevtORneppV]
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
    retrieveWf( wfs, w_cx, nevt, 4 );
    retrieveWf( wfs, w_cx, nevt, 6 );
    retrieveWf( wfs, w_cx, nevt, 7 );
    retrieveWf( wfs, w_cx, nevt, 80 );
    retrieveWf( wfs, w_cx, nevt, 83 );
    retrieveWf( wfs, w_cx, nevt, 84 );
    retrieveWf( wfs, w_cx, nevt, 86 );
    retrieveWf( wfs, w_cx, nevt, 98 );
    retrieveWf( wfs, w_cx, nevt, 102 );
    retrieveWf( wfs, w_cx, nevt, 118 );
    retrieveWf( wfs, w_cx, nevt, 128 );
    retrieveWf( wfs, w_cx, nevt, 129 );
    retrieveWf( wfs, w_cx, nevt, 144 );
    retrieveWf( wfs, w_cx, nevt, 150 );
    retrieveWf( wfs, w_cx, nevt, 174 );
    retrieveWf( wfs, w_cx, nevt, 176 );
    retrieveWf( wfs, w_cx, nevt, 179 );
    retrieveWf( wfs, w_cx, nevt, 181 );
    retrieveWf( wfs, w_cx, nevt, 196 );
    retrieveWf( wfs, w_cx, nevt, 202 );
    retrieveWf( wfs, w_cx, nevt, 204 );
    retrieveWf( wfs, w_cx, nevt, 239 );
    retrieveWf( wfs, w_cx, nevt, 270 );
    retrieveWf( wfs, w_cx, nevt, 271 );
    retrieveWf( wfs, w_cx, nevt, 279 );
    retrieveWf( wfs, w_cx, nevt, 283 );
    retrieveWf( wfs, w_cx, nevt, 286 );
    retrieveWf( wfs, w_cx, nevt, 289 );
    retrieveWf( wfs, w_cx, nevt, 290 );
    retrieveWf( wfs, w_cx, nevt, 292 );
    retrieveWf( wfs, w_cx, nevt, 293 );
    retrieveWf( wfs, w_cx, nevt, 294 );
    retrieveWf( wfs, w_cx, nevt, 295 );
    retrieveWf( wfs, w_cx, nevt, 296 );
    retrieveWf( wfs, w_cx, nevt, 297 );
    retrieveWf( wfs, w_cx, nevt, 298 );
    retrieveWf( wfs, w_cx, nevt, 306 );
    retrieveWf( wfs, w_cx, nevt, 307 );
    retrieveWf( wfs, w_cx, nevt, 308 );
    retrieveWf( wfs, w_cx, nevt, 309 );
    retrieveWf( wfs, w_cx, nevt, 310 );
    retrieveWf( wfs, w_cx, nevt, 311 );
    retrieveWf( wfs, w_cx, nevt, 435 );
    retrieveWf( wfs, w_cx, nevt, 493 );
    retrieveWf( wfs, w_cx, nevt, 494 );
    retrieveWf( wfs, w_cx, nevt, 495 );
    retrieveWf( wfs, w_cx, nevt, 498 );
    retrieveWf( wfs, w_cx, nevt, 500 );
    retrieveWf( wfs, w_cx, nevt, 516 );
    retrieveWf( wfs, w_cx, nevt, 518 );
    retrieveWf( wfs, w_cx, nevt, 536 );
    retrieveWf( wfs, w_cx, nevt, 544 );
    retrieveWf( wfs, w_cx, nevt, 545 );
    retrieveWf( wfs, w_cx, nevt, 555 );
    retrieveWf( wfs, w_cx, nevt, 576 );
    retrieveWf( wfs, w_cx, nevt, 578 );
    retrieveWf( wfs, w_cx, nevt, 583 );
    retrieveWf( wfs, w_cx, nevt, 584 );
    retrieveWf( wfs, w_cx, nevt, 594 );
    retrieveWf( wfs, w_cx, nevt, 595 );
    retrieveWf( wfs, w_cx, nevt, 596 );
    retrieveWf( wfs, w_cx, nevt, 662 );
    retrieveWf( wfs, w_cx, nevt, 670 );
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
    retrieveWf( wfs, w_cx, nevt, 724 );
    retrieveWf( wfs, w_cx, nevt, 725 );
    retrieveWf( wfs, w_cx, nevt, 727 );
    retrieveWf( wfs, w_cx, nevt, 728 );
    retrieveWf( wfs, w_cx, nevt, 734 );
    retrieveWf( wfs, w_cx, nevt, 735 );
    retrieveWf( wfs, w_cx, nevt, 743 );
    retrieveWf( wfs, w_cx, nevt, 748 );
#endif
#endif

    // *** DIAGRAM 10701 OF 15495 ***
    // Wavefunction(s) for diagram number 10701
    // (none)
    // Amplitude(s) for diagram number 10701
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[494], w_fp[84], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];

    // *** DIAGRAM 10702 OF 15495 ***
    // Wavefunction(s) for diagram number 10702
    // (none)
    // Amplitude(s) for diagram number 10702
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[2], w_fp[311], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[184] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[586] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[596] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[706] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[716] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10703 OF 15495 ***
    // Wavefunction(s) for diagram number 10703
    // (none)
    // Amplitude(s) for diagram number 10703
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[536], w_fp[128], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[586] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 10704 OF 15495 ***
    // Wavefunction(s) for diagram number 10704
    // (none)
    // Amplitude(s) for diagram number 10704
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[494], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[173] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[184] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[190] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[387] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[389] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[400] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[406] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10705 OF 15495 ***
    // Wavefunction(s) for diagram number 10705
    // (none)
    // Amplitude(s) for diagram number 10705
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[306], w_fp[2], w_fp[516], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[94] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[115] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[118] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[577] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[580] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[697] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[700] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10706 OF 15495 ***
    // Wavefunction(s) for diagram number 10706
    // (none)
    // Amplitude(s) for diagram number 10706
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[516], w_fp[279], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] -= amp_sv[0];
    jamp_sv[94] += amp_sv[0];
    jamp_sv[115] += amp_sv[0];
    jamp_sv[118] -= amp_sv[0];
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[184] -= amp_sv[0];
    jamp_sv[190] += amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];
    jamp_sv[400] += amp_sv[0];
    jamp_sv[406] -= amp_sv[0];
    jamp_sv[577] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 10707 OF 15495 ***
    // Wavefunction(s) for diagram number 10707
    // (none)
    // Amplitude(s) for diagram number 10707
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[129], w_fp[494], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[171] += amp_sv[0];
    jamp_sv[173] -= amp_sv[0];
    jamp_sv[387] -= amp_sv[0];
    jamp_sv[389] += amp_sv[0];

    // *** DIAGRAM 10708 OF 15495 ***
    // Wavefunction(s) for diagram number 10708
    // (none)
    // Amplitude(s) for diagram number 10708
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[306], w_fp[128], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[577] += amp_sv[0];
    jamp_sv[580] -= amp_sv[0];
    jamp_sv[697] -= amp_sv[0];
    jamp_sv[700] += amp_sv[0];

    // *** DIAGRAM 10709 OF 15495 ***
    // Wavefunction(s) for diagram number 10709
    // (none)
    // Amplitude(s) for diagram number 10709
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[311], w_fp[239], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];

    // *** DIAGRAM 10710 OF 15495 ***
    // Wavefunction(s) for diagram number 10710
    // (none)
    // Amplitude(s) for diagram number 10710
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[286], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] -= amp_sv[0];
    jamp_sv[11] += amp_sv[0];
    jamp_sv[51] += amp_sv[0];
    jamp_sv[53] -= amp_sv[0];
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[586] -= amp_sv[0];
    jamp_sv[596] += amp_sv[0];
    jamp_sv[706] += amp_sv[0];
    jamp_sv[716] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[727], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[91] += amp_sv[0];
    jamp_sv[94] -= amp_sv[0];
    jamp_sv[115] -= amp_sv[0];
    jamp_sv[118] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[184] += amp_sv[0];
    jamp_sv[190] -= amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[400] -= amp_sv[0];
    jamp_sv[406] += amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[196], w_fp[2], w_fp[728], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[9] += amp_sv[0];
    jamp_sv[11] -= amp_sv[0];
    jamp_sv[51] -= amp_sv[0];
    jamp_sv[53] += amp_sv[0];
    jamp_sv[171] -= amp_sv[0];
    jamp_sv[173] += amp_sv[0];
    jamp_sv[387] += amp_sv[0];
    jamp_sv[389] -= amp_sv[0];
    jamp_sv[577] -= amp_sv[0];
    jamp_sv[580] += amp_sv[0];
    jamp_sv[586] += amp_sv[0];
    jamp_sv[596] -= amp_sv[0];
    jamp_sv[697] += amp_sv[0];
    jamp_sv[700] -= amp_sv[0];
    jamp_sv[706] -= amp_sv[0];
    jamp_sv[716] += amp_sv[0];

    // *** DIAGRAM 10711 OF 15495 ***
    // Wavefunction(s) for diagram number 10711
    // (none)
    // Amplitude(s) for diagram number 10711
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[4], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 10712 OF 15495 ***
    // Wavefunction(s) for diagram number 10712
    // (none)
    // Amplitude(s) for diagram number 10712
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[7], w_fp[706], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 10713 OF 15495 ***
    // Wavefunction(s) for diagram number 10713
    // (none)
    // Amplitude(s) for diagram number 10713
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[4], w_fp[707], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 10714 OF 15495 ***
    // Wavefunction(s) for diagram number 10714
    // (none)
    // Amplitude(s) for diagram number 10714
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[500], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];

    // *** DIAGRAM 10715 OF 15495 ***
    // Wavefunction(s) for diagram number 10715
    // (none)
    // Amplitude(s) for diagram number 10715
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[707], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10716 OF 15495 ***
    // Wavefunction(s) for diagram number 10716
    // (none)
    // Amplitude(s) for diagram number 10716
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[500], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];

    // *** DIAGRAM 10717 OF 15495 ***
    // Wavefunction(s) for diagram number 10717
    // (none)
    // Amplitude(s) for diagram number 10717
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10718 OF 15495 ***
    // Wavefunction(s) for diagram number 10718
    // (none)
    // Amplitude(s) for diagram number 10718
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[595], w_fp[494], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10719 OF 15495 ***
    // Wavefunction(s) for diagram number 10719
    // (none)
    // Amplitude(s) for diagram number 10719
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[578], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10720 OF 15495 ***
    // Wavefunction(s) for diagram number 10720
    // (none)
    // Amplitude(s) for diagram number 10720
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[7], w_fp[724], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10721 OF 15495 ***
    // Wavefunction(s) for diagram number 10721
    // (none)
    // Amplitude(s) for diagram number 10721
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[578], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += amp_sv[0];
    jamp_sv[274] -= amp_sv[0];
    jamp_sv[298] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];

    // *** DIAGRAM 10722 OF 15495 ***
    // Wavefunction(s) for diagram number 10722
    // (none)
    // Amplitude(s) for diagram number 10722
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[292], w_fp[4], w_fp[724], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10723 OF 15495 ***
    // Wavefunction(s) for diagram number 10723
    // (none)
    // Amplitude(s) for diagram number 10723
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[595], w_fp[2], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[638] -= amp_sv[0];
    jamp_sv[680] += amp_sv[0];

    // *** DIAGRAM 10724 OF 15495 ***
    // Wavefunction(s) for diagram number 10724
    // (none)
    // Amplitude(s) for diagram number 10724
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[296], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[297], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[298], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[274] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[298] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[638] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[680] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10725 OF 15495 ***
    // Wavefunction(s) for diagram number 10725
    // (none)
    // Amplitude(s) for diagram number 10725
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[498], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10726 OF 15495 ***
    // Wavefunction(s) for diagram number 10726
    // (none)
    // Amplitude(s) for diagram number 10726
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[494], w_fp[7], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10727 OF 15495 ***
    // Wavefunction(s) for diagram number 10727
    // (none)
    // Amplitude(s) for diagram number 10727
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[498], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10728 OF 15495 ***
    // Wavefunction(s) for diagram number 10728
    // (none)
    // Amplitude(s) for diagram number 10728
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[555], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10729 OF 15495 ***
    // Wavefunction(s) for diagram number 10729
    // (none)
    // Amplitude(s) for diagram number 10729
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[289], w_fp[150], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[289], w_fp[150], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[289], w_fp[150], w_fp[7], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[175] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[391] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[606] -= amp_sv[0];
    jamp_sv[608] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];

    // *** DIAGRAM 10730 OF 15495 ***
    // Wavefunction(s) for diagram number 10730
    // (none)
    // Amplitude(s) for diagram number 10730
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[7], w_fp[743], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 10731 OF 15495 ***
    // Wavefunction(s) for diagram number 10731
    // (none)
    // Amplitude(s) for diagram number 10731
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[7], w_fp[725], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 10732 OF 15495 ***
    // Wavefunction(s) for diagram number 10732
    // (none)
    // Amplitude(s) for diagram number 10732
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[743], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10733 OF 15495 ***
    // Wavefunction(s) for diagram number 10733
    // (none)
    // Amplitude(s) for diagram number 10733
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[555], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[175] += amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[391] -= amp_sv[0];

    // *** DIAGRAM 10734 OF 15495 ***
    // Wavefunction(s) for diagram number 10734
    // (none)
    // Amplitude(s) for diagram number 10734
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[292], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[292], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[292], w_fp[150], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[186] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[258] -= amp_sv[0];
    jamp_sv[260] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 10735 OF 15495 ***
    // Wavefunction(s) for diagram number 10735
    // (none)
    // Amplitude(s) for diagram number 10735
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[150], w_fp[4], w_fp[493], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];

    // *** DIAGRAM 10736 OF 15495 ***
    // Wavefunction(s) for diagram number 10736
    // (none)
    // Amplitude(s) for diagram number 10736
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[292], w_fp[4], w_fp[725], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 10737 OF 15495 ***
    // Wavefunction(s) for diagram number 10737
    // (none)
    // Amplitude(s) for diagram number 10737
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[493], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10738 OF 15495 ***
    // Wavefunction(s) for diagram number 10738
    // (none)
    // Amplitude(s) for diagram number 10738
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[545], w_fp[2], w_fp[292], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[186] += amp_sv[0];
    jamp_sv[402] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 10739 OF 15495 ***
    // Wavefunction(s) for diagram number 10739
    // (none)
    // Amplitude(s) for diagram number 10739
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[544], w_fp[150], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[271] += amp_sv[0];
    jamp_sv[295] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[602] -= amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[612] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[711], w_fp[150], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[175] -= amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[391] += amp_sv[0];
    jamp_sv[606] += amp_sv[0];
    jamp_sv[608] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[678] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[710], w_fp[150], w_fp[7], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 10740 OF 15495 ***
    // Wavefunction(s) for diagram number 10740
    // (none)
    // Amplitude(s) for diagram number 10740
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[25] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[31] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[175] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[391] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[176], w_fp[2], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[241] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[247] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[271] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[295] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10741 OF 15495 ***
    // Wavefunction(s) for diagram number 10741
    // (none)
    // Amplitude(s) for diagram number 10741
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[714], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[241] -= amp_sv[0];
    jamp_sv[247] += amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[636] += amp_sv[0];
    jamp_sv[678] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[713], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[186] -= amp_sv[0];
    jamp_sv[258] += amp_sv[0];
    jamp_sv[260] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[402] += amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[712], w_fp[150], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[241] += amp_sv[0];
    jamp_sv[247] -= amp_sv[0];
    jamp_sv[271] -= amp_sv[0];
    jamp_sv[295] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[602] += amp_sv[0];
    jamp_sv[612] -= amp_sv[0];
    jamp_sv[636] -= amp_sv[0];
    jamp_sv[678] += amp_sv[0];

    // *** DIAGRAM 10742 OF 15495 ***
    // Wavefunction(s) for diagram number 10742
    // (none)
    // Amplitude(s) for diagram number 10742
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[714], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[713], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[98] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[108] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[186] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[402] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[202], w_fp[2], w_fp[712], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[602] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[612] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[636] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[678] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10743 OF 15495 ***
    // Wavefunction(s) for diagram number 10743
    // (none)
    // Amplitude(s) for diagram number 10743
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[296], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[297], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[10] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[298], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[25] -= amp_sv[0];
    jamp_sv[31] += amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[98] -= amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[108] += amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[274] += amp_sv[0];
    jamp_sv[298] -= amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[638] += amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];
    jamp_sv[680] -= amp_sv[0];

    // *** DIAGRAM 10744 OF 15495 ***
    // Wavefunction(s) for diagram number 10744
    // (none)
    // Amplitude(s) for diagram number 10744
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[150], w_fp[102], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];

    // *** DIAGRAM 10745 OF 15495 ***
    // Wavefunction(s) for diagram number 10745
    // (none)
    // Amplitude(s) for diagram number 10745
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[98], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10746 OF 15495 ***
    // Wavefunction(s) for diagram number 10746
    // (none)
    // Amplitude(s) for diagram number 10746
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[83], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[10] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[49] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[52] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[169] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10747 OF 15495 ***
    // Wavefunction(s) for diagram number 10747
    // (none)
    // Amplitude(s) for diagram number 10747
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[494], w_fp[102], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];

    // *** DIAGRAM 10748 OF 15495 ***
    // Wavefunction(s) for diagram number 10748
    // (none)
    // Amplitude(s) for diagram number 10748
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[2], w_fp[310], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[178] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[344] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[350] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[656] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[662] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10749 OF 15495 ***
    // Wavefunction(s) for diagram number 10749
    // (none)
    // Amplitude(s) for diagram number 10749
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[594], w_fp[98], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];

    // *** DIAGRAM 10750 OF 15495 ***
    // Wavefunction(s) for diagram number 10750
    // (none)
    // Amplitude(s) for diagram number 10750
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[494], w_fp[518], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[172] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[178] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[188] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[385] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[388] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[394] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[404] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10751 OF 15495 ***
    // Wavefunction(s) for diagram number 10751
    // (none)
    // Amplitude(s) for diagram number 10751
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[307], w_fp[2], w_fp[518], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[44] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[102] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[104] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[336] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[338] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[648] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[650] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10752 OF 15495 ***
    // Wavefunction(s) for diagram number 10752
    // (none)
    // Amplitude(s) for diagram number 10752
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[518], w_fp[279], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] -= amp_sv[0];
    jamp_sv[44] += amp_sv[0];
    jamp_sv[102] += amp_sv[0];
    jamp_sv[104] -= amp_sv[0];
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[178] -= amp_sv[0];
    jamp_sv[188] += amp_sv[0];
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];
    jamp_sv[394] += amp_sv[0];
    jamp_sv[404] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];

    // *** DIAGRAM 10753 OF 15495 ***
    // Wavefunction(s) for diagram number 10753
    // (none)
    // Amplitude(s) for diagram number 10753
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[83], w_fp[494], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[169] += amp_sv[0];
    jamp_sv[172] -= amp_sv[0];
    jamp_sv[385] -= amp_sv[0];
    jamp_sv[388] += amp_sv[0];

    // *** DIAGRAM 10754 OF 15495 ***
    // Wavefunction(s) for diagram number 10754
    // (none)
    // Amplitude(s) for diagram number 10754
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[307], w_fp[98], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[336] += amp_sv[0];
    jamp_sv[338] -= amp_sv[0];
    jamp_sv[648] -= amp_sv[0];
    jamp_sv[650] += amp_sv[0];

    // *** DIAGRAM 10755 OF 15495 ***
    // Wavefunction(s) for diagram number 10755
    // (none)
    // Amplitude(s) for diagram number 10755
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[310], w_fp[150], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];

    // *** DIAGRAM 10756 OF 15495 ***
    // Wavefunction(s) for diagram number 10756
    // (none)
    // Amplitude(s) for diagram number 10756
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[283], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] -= amp_sv[0];
    jamp_sv[10] += amp_sv[0];
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[49] += amp_sv[0];
    jamp_sv[52] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[344] -= amp_sv[0];
    jamp_sv[350] += amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[656] += amp_sv[0];
    jamp_sv[662] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[271], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[42] += amp_sv[0];
    jamp_sv[44] -= amp_sv[0];
    jamp_sv[102] -= amp_sv[0];
    jamp_sv[104] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[178] += amp_sv[0];
    jamp_sv[188] -= amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[394] -= amp_sv[0];
    jamp_sv[404] += amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[174], w_fp[2], w_fp[270], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[7] += amp_sv[0];
    jamp_sv[10] -= amp_sv[0];
    jamp_sv[49] -= amp_sv[0];
    jamp_sv[52] += amp_sv[0];
    jamp_sv[169] -= amp_sv[0];
    jamp_sv[172] += amp_sv[0];
    jamp_sv[336] -= amp_sv[0];
    jamp_sv[338] += amp_sv[0];
    jamp_sv[344] += amp_sv[0];
    jamp_sv[350] -= amp_sv[0];
    jamp_sv[385] += amp_sv[0];
    jamp_sv[388] -= amp_sv[0];
    jamp_sv[648] += amp_sv[0];
    jamp_sv[650] -= amp_sv[0];
    jamp_sv[656] -= amp_sv[0];
    jamp_sv[662] += amp_sv[0];

    // *** DIAGRAM 10757 OF 15495 ***
    // Wavefunction(s) for diagram number 10757
    // (none)
    // Amplitude(s) for diagram number 10757
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[144], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
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
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[144], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[144], w_fp[4], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 10758 OF 15495 ***
    // Wavefunction(s) for diagram number 10758
    // (none)
    // Amplitude(s) for diagram number 10758
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[6], w_fp[706], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[528] += amp_sv[0];
    jamp_sv[530] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 10759 OF 15495 ***
    // Wavefunction(s) for diagram number 10759
    // (none)
    // Amplitude(s) for diagram number 10759
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[4], w_fp[708], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 10760 OF 15495 ***
    // Wavefunction(s) for diagram number 10760
    // (none)
    // Amplitude(s) for diagram number 10760
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[500], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];

    // *** DIAGRAM 10761 OF 15495 ***
    // Wavefunction(s) for diagram number 10761
    // (none)
    // Amplitude(s) for diagram number 10761
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[2], w_fp[708], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10762 OF 15495 ***
    // Wavefunction(s) for diagram number 10762
    // (none)
    // Amplitude(s) for diagram number 10762
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[500], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[168] -= amp_sv[0];
    jamp_sv[384] += amp_sv[0];

    // *** DIAGRAM 10763 OF 15495 ***
    // Wavefunction(s) for diagram number 10763
    // (none)
    // Amplitude(s) for diagram number 10763
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[706], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10764 OF 15495 ***
    // Wavefunction(s) for diagram number 10764
    // (none)
    // Amplitude(s) for diagram number 10764
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[494], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10765 OF 15495 ***
    // Wavefunction(s) for diagram number 10765
    // (none)
    // Amplitude(s) for diagram number 10765
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10766 OF 15495 ***
    // Wavefunction(s) for diagram number 10766
    // (none)
    // Amplitude(s) for diagram number 10766
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[6], w_fp[735], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10767 OF 15495 ***
    // Wavefunction(s) for diagram number 10767
    // (none)
    // Amplitude(s) for diagram number 10767
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[495], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += amp_sv[0];
    jamp_sv[272] -= amp_sv[0];
    jamp_sv[296] += amp_sv[0];
    jamp_sv[392] -= amp_sv[0];

    // *** DIAGRAM 10768 OF 15495 ***
    // Wavefunction(s) for diagram number 10768
    // (none)
    // Amplitude(s) for diagram number 10768
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[290], w_fp[4], w_fp[735], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10769 OF 15495 ***
    // Wavefunction(s) for diagram number 10769
    // (none)
    // Amplitude(s) for diagram number 10769
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[435], w_fp[2], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[518] -= amp_sv[0];
    jamp_sv[560] += amp_sv[0];

    // *** DIAGRAM 10770 OF 15495 ***
    // Wavefunction(s) for diagram number 10770
    // (none)
    // Amplitude(s) for diagram number 10770
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[293], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[294], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[295], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[272] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[296] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[518] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[560] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10771 OF 15495 ***
    // Wavefunction(s) for diagram number 10771
    // (none)
    // Amplitude(s) for diagram number 10771
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[498], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10772 OF 15495 ***
    // Wavefunction(s) for diagram number 10772
    // (none)
    // Amplitude(s) for diagram number 10772
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[494], w_fp[6], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10773 OF 15495 ***
    // Wavefunction(s) for diagram number 10773
    // (none)
    // Amplitude(s) for diagram number 10773
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[498], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10774 OF 15495 ***
    // Wavefunction(s) for diagram number 10774
    // (none)
    // Amplitude(s) for diagram number 10774
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[494], w_fp[4], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10775 OF 15495 ***
    // Wavefunction(s) for diagram number 10775
    // (none)
    // Amplitude(s) for diagram number 10775
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[289], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[289], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[289], w_fp[144], w_fp[6], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[174] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[390] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[486] -= amp_sv[0];
    jamp_sv[488] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];

    // *** DIAGRAM 10776 OF 15495 ***
    // Wavefunction(s) for diagram number 10776
    // (none)
    // Amplitude(s) for diagram number 10776
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[6], w_fp[743], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 10777 OF 15495 ***
    // Wavefunction(s) for diagram number 10777
    // (none)
    // Amplitude(s) for diagram number 10777
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[289], w_fp[6], w_fp[734], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 10778 OF 15495 ***
    // Wavefunction(s) for diagram number 10778
    // (none)
    // Amplitude(s) for diagram number 10778
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[743], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10779 OF 15495 ***
    // Wavefunction(s) for diagram number 10779
    // (none)
    // Amplitude(s) for diagram number 10779
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[583], w_fp[2], w_fp[289], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[174] += amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[390] -= amp_sv[0];

    // *** DIAGRAM 10780 OF 15495 ***
    // Wavefunction(s) for diagram number 10780
    // (none)
    // Amplitude(s) for diagram number 10780
    VVVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[290], w_fp[144], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    VVVV3_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[290], w_fp[144], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    VVVV4_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[290], w_fp[144], w_fp[4], COUPs[2], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[180] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[252] -= amp_sv[0];
    jamp_sv[254] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 10781 OF 15495 ***
    // Wavefunction(s) for diagram number 10781
    // (none)
    // Amplitude(s) for diagram number 10781
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[144], w_fp[4], w_fp[748], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];

    // *** DIAGRAM 10782 OF 15495 ***
    // Wavefunction(s) for diagram number 10782
    // (none)
    // Amplitude(s) for diagram number 10782
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[290], w_fp[4], w_fp[734], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 10783 OF 15495 ***
    // Wavefunction(s) for diagram number 10783
    // (none)
    // Amplitude(s) for diagram number 10783
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[2], w_fp[748], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10784 OF 15495 ***
    // Wavefunction(s) for diagram number 10784
    // (none)
    // Amplitude(s) for diagram number 10784
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[584], w_fp[2], w_fp[290], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[180] += amp_sv[0];
    jamp_sv[396] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];

    // *** DIAGRAM 10785 OF 15495 ***
    // Wavefunction(s) for diagram number 10785
    // (none)
    // Amplitude(s) for diagram number 10785
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[544], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[270] += amp_sv[0];
    jamp_sv[294] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[482] -= amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[492] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[711], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[174] -= amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[390] += amp_sv[0];
    jamp_sv[486] += amp_sv[0];
    jamp_sv[488] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[558] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[710], w_fp[144], w_fp[6], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= amp_sv[0];
    jamp_sv[48] += amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];
    jamp_sv[536] += amp_sv[0];
    jamp_sv[542] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];

    // *** DIAGRAM 10786 OF 15495 ***
    // Wavefunction(s) for diagram number 10786
    // (none)
    // Amplitude(s) for diagram number 10786
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[544], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[711], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[24] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[30] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[174] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[390] -= cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[181], w_fp[2], w_fp[710], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[240] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[246] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[270] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[294] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10787 OF 15495 ***
    // Wavefunction(s) for diagram number 10787
    // (none)
    // Amplitude(s) for diagram number 10787
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[709], w_fp[144], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[240] -= amp_sv[0];
    jamp_sv[246] += amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[516] += amp_sv[0];
    jamp_sv[558] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[703], w_fp[144], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[180] -= amp_sv[0];
    jamp_sv[252] += amp_sv[0];
    jamp_sv[254] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[396] += amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[670], w_fp[144], w_fp[4], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= amp_sv[0];
    jamp_sv[50] += amp_sv[0];
    jamp_sv[170] += amp_sv[0];
    jamp_sv[240] += amp_sv[0];
    jamp_sv[246] -= amp_sv[0];
    jamp_sv[270] -= amp_sv[0];
    jamp_sv[294] += amp_sv[0];
    jamp_sv[312] -= amp_sv[0];
    jamp_sv[314] += amp_sv[0];
    jamp_sv[320] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[386] -= amp_sv[0];
    jamp_sv[482] += amp_sv[0];
    jamp_sv[492] -= amp_sv[0];
    jamp_sv[516] -= amp_sv[0];
    jamp_sv[558] += amp_sv[0];

    // *** DIAGRAM 10788 OF 15495 ***
    // Wavefunction(s) for diagram number 10788
    // (none)
    // Amplitude(s) for diagram number 10788
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[2], w_fp[709], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[2], w_fp[703], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[74] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[84] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[180] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[396] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[204], w_fp[2], w_fp[670], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[482] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[492] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[516] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[558] -= cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10789 OF 15495 ***
    // Wavefunction(s) for diagram number 10789
    // (none)
    // Amplitude(s) for diagram number 10789
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[293], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
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
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[294], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[8] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[36] += amp_sv[0];
    jamp_sv[38] -= amp_sv[0];
    jamp_sv[50] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[320] -= amp_sv[0];
    jamp_sv[326] += amp_sv[0];
    jamp_sv[398] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[0], w_fp[295], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] += amp_sv[0];
    jamp_sv[24] -= amp_sv[0];
    jamp_sv[30] += amp_sv[0];
    jamp_sv[48] -= amp_sv[0];
    jamp_sv[74] -= amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[84] += amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[272] += amp_sv[0];
    jamp_sv[296] -= amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[518] += amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];
    jamp_sv[560] -= amp_sv[0];

    // *** DIAGRAM 10790 OF 15495 ***
    // Wavefunction(s) for diagram number 10790
    // (none)
    // Amplitude(s) for diagram number 10790
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[662], w_fp[144], w_fp[86], COUPs[0], 1.0, &amp_fp[0] );
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

    // *** DIAGRAM 10791 OF 15495 ***
    // Wavefunction(s) for diagram number 10791
    // (none)
    // Amplitude(s) for diagram number 10791
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[118], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10792 OF 15495 ***
    // Wavefunction(s) for diagram number 10792
    // (none)
    // Amplitude(s) for diagram number 10792
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[80], w_fp[2], w_fp[662], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[6] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[8] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[48] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[50] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[168] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10793 OF 15495 ***
    // Wavefunction(s) for diagram number 10793
    // (none)
    // Amplitude(s) for diagram number 10793
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[494], w_fp[86], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] += amp_sv[0];
    jamp_sv[182] -= amp_sv[0];
    jamp_sv[392] -= amp_sv[0];
    jamp_sv[398] += amp_sv[0];

    // *** DIAGRAM 10794 OF 15495 ***
    // Wavefunction(s) for diagram number 10794
    // (none)
    // Amplitude(s) for diagram number 10794
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[2], w_fp[309], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[176] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[320] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[326] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[536] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[542] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10795 OF 15495 ***
    // Wavefunction(s) for diagram number 10795
    // (none)
    // Amplitude(s) for diagram number 10795
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[596], w_fp[118], w_fp[279], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[320] += amp_sv[0];
    jamp_sv[326] -= amp_sv[0];
    jamp_sv[536] -= amp_sv[0];
    jamp_sv[542] += amp_sv[0];

    // *** DIAGRAM 10796 OF 15495 ***
    // Wavefunction(s) for diagram number 10796
    // (none)
    // Amplitude(s) for diagram number 10796
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[179], w_fp[494], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[170] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[176] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[182] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[384] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[386] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[392] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[398] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10797 OF 15495 ***
    // Wavefunction(s) for diagram number 10797
    // (none)
    // Amplitude(s) for diagram number 10797
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[308], w_fp[2], w_fp[576], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[38] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[78] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[80] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[312] += cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[314] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[528] -= cxtype( 0, 1 ) * amp_sv[0];
    jamp_sv[530] += cxtype( 0, 1 ) * amp_sv[0];

    // *** DIAGRAM 10798 OF 15495 ***
    // Wavefunction(s) for diagram number 10798
    // (none)
    // Amplitude(s) for diagram number 10798
    VVV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[576], w_fp[279], w_fp[144], COUPs[0], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[36] -= amp_sv[0];
    jamp_sv[38] += amp_sv[0];
    jamp_sv[78] += amp_sv[0];
    jamp_sv[80] -= amp_sv[0];
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[176] -= amp_sv[0];
    jamp_sv[182] += amp_sv[0];
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];
    jamp_sv[392] += amp_sv[0];
    jamp_sv[398] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];

    // *** DIAGRAM 10799 OF 15495 ***
    // Wavefunction(s) for diagram number 10799
    // (none)
    // Amplitude(s) for diagram number 10799
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[80], w_fp[494], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[168] += amp_sv[0];
    jamp_sv[170] -= amp_sv[0];
    jamp_sv[384] -= amp_sv[0];
    jamp_sv[386] += amp_sv[0];

    // *** DIAGRAM 10800 OF 15495 ***
    // Wavefunction(s) for diagram number 10800
    // (none)
    // Amplitude(s) for diagram number 10800
    FFV1_0<W_ACCESS, A_ACCESS, CD_ACCESS>( w_fp[308], w_fp[118], w_fp[0], COUPs[1], 1.0, &amp_fp[0] );
#ifdef MGONGPU_SUPPORTS_MULTICHANNEL
    // Here the code base generated with multichannel support updates numerators_sv and denominators_sv (#473)
#endif
    jamp_sv[312] += amp_sv[0];
    jamp_sv[314] -= amp_sv[0];
    jamp_sv[528] -= amp_sv[0];
    jamp_sv[530] += amp_sv[0];

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
